import json
from datetime import datetime

import pandas as pd

from covigator.configuration import Configuration
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorErrorProcessingCoverageResults, CovigatorExcludedSampleBadQualityReads, \
    CovigatorExcludedSampleNarrowCoverage
from covigator.misc import backoff_retrier
from covigator.database.model import JobStatus, JobEna, Sample, DataSource
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.cooccurrence_matrix import CooccurrenceMatrix
from covigator.pipeline.downloader import Downloader
from covigator.pipeline.ena_pipeline import Pipeline
from covigator.pipeline.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 10


class EnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(EnaProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA,
            function=EnaProcessor.run_all)

    @staticmethod
    def run_all(job: JobEna, queries: Queries, config: Configuration):
        EnaProcessor.download(job=job, queries=queries, config=config)
        EnaProcessor.run_pipeline(job=job, queries=queries, config=config)
        EnaProcessor.load(job=job, queries=queries, config=config)
        EnaProcessor.compute_cooccurrence(job=job, queries=queries, config=config)

    @staticmethod
    def download(job: JobEna, queries: Queries, config: Configuration):
        # ensures that the download is done with retries, even after MD5 check sum failure
        downloader = Downloader(config=config)
        download_with_retries = backoff_retrier.wrapper(downloader.download, NUMBER_RETRIES_DOWNLOADER)
        sample_ena = queries.find_sample_by_accession(job.run_accession, source=DataSource.ENA)
        paths = download_with_retries(sample_ena=sample_ena)
        job.fastq_path = paths
        job.downloaded_at = datetime.now()

    @staticmethod
    def run_pipeline(job: JobEna, queries: Queries, config: Configuration):
        fastq1, fastq2 = job.get_fastq1_and_fastq2()
        vcf_path, qc_path, vertical_coverage_path, horizontal_coverage_path = Pipeline(config=config)\
            .run(run_accession=job.run_accession, fastq1=fastq1, fastq2=fastq2)
        job.analysed_at = datetime.now()
        job.vcf_path = vcf_path
        job.qc_path = qc_path
        job.qc = json.load(open(qc_path))
        job.horizontal_coverage_path = horizontal_coverage_path
        job.vertical_coverage_path = vertical_coverage_path
        EnaProcessor.load_coverage_results(horizontal_coverage_path, job)

    @staticmethod
    def load_coverage_results(horizontal_coverage_path, job):
        try:
            data = pd.read_csv(horizontal_coverage_path, sep="\t")
            job.mean_depth = float(data.meandepth.loc[0])
            job.mean_base_quality = float(data.meanbaseq.loc[0])
            job.mean_mapping_quality = float(data.meanmapq.loc[0])
            job.num_reads = int(data.numreads.loc[0])
            job.num_reads = int(data.numreads.loc[0])
            job.covered_bases = int(data.covbases.loc[0])
            job.coverage = float(data.coverage.loc[0])
        except Exception as e:
            raise CovigatorErrorProcessingCoverageResults(e)

    @staticmethod
    def load(job: JobEna, queries: Queries, config: Configuration):
        if job.mean_mapping_quality < 10 or job.mean_base_quality < 10:
            raise CovigatorExcludedSampleBadQualityReads("Mean MQ: {}; mean BCQ: {}".format(
                job.mean_mapping_quality, job.mean_base_quality))
        if job.coverage < 20.0:
            raise CovigatorExcludedSampleNarrowCoverage("Horizontal coverage {} %".format(job.coverage))
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.loaded_at = datetime.now()

    @staticmethod
    def compute_cooccurrence(job: JobEna, queries: Queries, config: Configuration):
        CooccurrenceMatrix().compute(
            sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.cooccurrence_at = datetime.now()
