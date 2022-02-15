import json
import numpy as np
from datetime import datetime

import pandas as pd

from covigator.configuration import Configuration
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorErrorProcessingCoverageResults, CovigatorExcludedSampleBadQualityReads, \
    CovigatorExcludedSampleNarrowCoverage, CovigatorErrorProcessingPangolinResults, \
    CovigatorErrorProcessingDeduplicationResults
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

    def __init__(self, database: Database, dask_client: Client, config: Configuration, wait_time=60):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config, wait_time=wait_time)

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
        pipeline_result = Pipeline(config=config)\
            .run(run_accession=job.run_accession, fastq1=fastq1, fastq2=fastq2)
        job.analysed_at = datetime.now()

        # stores the paths to all files output by pipeline
        job.lofreq_vcf_path = pipeline_result.lofreq_vcf
        job.ivar_vcf_path = pipeline_result.ivar_vcf
        job.gatk_vcf_path = pipeline_result.gatk_vcf
        job.bcftools_vcf_path = pipeline_result.bcftools_vcf
        job.lofreq_pangolin_path = pipeline_result.lofreq_pangolin
        job.ivar_pangolin_path = pipeline_result.ivar_pangolin
        job.gatk_pangolin_path = pipeline_result.gatk_pangolin
        job.bcftools_pangolin_path = pipeline_result.bcftools_pangolin
        job.fastp_path = pipeline_result.fastp_qc
        job.horizontal_coverage_path = pipeline_result.horizontal_coverage
        job.vertical_coverage_path = pipeline_result.vertical_coverage
        job.deduplication_metrics_path = pipeline_result.deduplication_metrics

        # load FAST JSON into the DB
        job.qc = json.load(open(pipeline_result.fastp_qc))
        # load horizontal coverage values in the database
        EnaProcessor.load_coverage_results(job)
        # load deduplication metrics
        EnaProcessor.load_deduplication_metrics(job)
        # load pangolin results
        EnaProcessor.load_lofreq_pangolin(job)

    @staticmethod
    def load_deduplication_metrics(job: JobEna):
        try:
            data = pd.read_csv(job.deduplication_metrics_path,
                               sep="\t",
                               skiprows=6,
                               dtype={
                                   'PERCENT_DUPLICATION': float,
                                   'UNPAIRED_READS_EXAMINED': int,
                                   'READ_PAIRS_EXAMINED': int,
                                   'SECONDARY_OR_SUPPLEMENTARY_RDS': int,
                                   'UNMAPPED_READS': int,
                                   'UNPAIRED_READ_DUPLICATES': int,
                                   'READ_PAIR_DUPLICATES': int,
                                   'READ_PAIR_OPTICAL_DUPLICATES': int
                               })
            data.fillna(0, inplace=True)
            job.percent_duplication = data.PERCENT_DUPLICATION.loc[0]
            job.unpaired_reads_examined = data.UNPAIRED_READS_EXAMINED.loc[0]
            job.read_pairs_examined = data.READ_PAIRS_EXAMINED.loc[0]
            job.secondary_or_supplementary_reads = data.SECONDARY_OR_SUPPLEMENTARY_RDS.loc[0]
            job.unmapped_reads = data.UNMAPPED_READS.loc[0]
            job.unpaired_read_duplicates = data.UNPAIRED_READ_DUPLICATES.loc[0]
            job.read_pair_duplicates = data.READ_PAIR_DUPLICATES.loc[0]
            job.read_pair_optical_duplicates = data.READ_PAIR_OPTICAL_DUPLICATES.loc[0]
        except Exception as e:
            raise CovigatorErrorProcessingDeduplicationResults(e)

    @staticmethod
    def load_lofreq_pangolin(job: JobEna):
        try:
            data = pd.read_csv(job.lofreq_pangolin_path,
                               na_values=None,
                               dtype={
                                   'lineage': str,
                                   'conflict': float,
                                   'ambiguity_score': float,
                                   'scorpio_call': str,
                                   'scorpio_support': float,
                                   'scorpio_conflict': float,
                                   'version': str,
                                   'pangolin_version': str,
                                   'pangoLEARN_version': str,
                                   'pango_version': str,
                                   'status': str,
                                   'note': str
                               })
            data.fillna(value="", inplace=True)
            job.pangolin_lineage = data.lineage.loc[0]
            job.pangolin_conflict = data.conflict.loc[0]
            job.pangolin_ambiguity_score = data.ambiguity_score.loc[0]
            job.pangolin_scorpio_call = data.scorpio_call.loc[0]
            job.pangolin_scorpio_support = data.scorpio_support.loc[0]
            job.pangolin_scorpio_conflict = data.scorpio_conflict.loc[0]
            job.pangolin_version = data.version.loc[0]
            job.pangolin_pangolin_version = data.pangolin_version.loc[0]
            job.pangolin_pangoLEARN_version = data.pangoLEARN_version.loc[0]
            job.pangolin_pango_version = data.pango_version.loc[0]
            job.pangolin_status = data.status.loc[0]
            job.pangolin_note = data.note.loc[0]
        except Exception as e:
            raise CovigatorErrorProcessingPangolinResults(e)

    @staticmethod
    def load_coverage_results(job: JobEna):
        try:
            data = pd.read_csv(job.horizontal_coverage_path,
                               sep="\t",
                               dtype={
                                   'numreads': int,
                                   'covbases': int,
                                   'meandepth': float,
                                   'meanbaseq': float,
                                   'meanmapq': float,
                                   'coverage': float,
                               }
                               )
            data.fillna(0, inplace=True)
            job.mean_depth = float(data.meandepth.loc[0])
            job.mean_base_quality = float(data.meanbaseq.loc[0])
            job.mean_mapping_quality = float(data.meanmapq.loc[0])
            job.num_reads = int(data.numreads.loc[0])
            job.covered_bases = int(data.covbases.loc[0])
            job.coverage = float(data.coverage.loc[0])
        except Exception as e:
            raise CovigatorErrorProcessingCoverageResults(e)

    @staticmethod
    def load(job: JobEna, queries: Queries, config: Configuration):
        if job.mean_mapping_quality < config.mean_mq_thr or job.mean_base_quality < config.mean_bq_thr:
            raise CovigatorExcludedSampleBadQualityReads("Mean MQ: {}; mean BCQ: {}".format(
                job.mean_mapping_quality, job.mean_base_quality))
        if job.coverage < config.horizontal_coverage_thr:
            raise CovigatorExcludedSampleNarrowCoverage("Horizontal coverage {} %".format(job.coverage))
        VcfLoader().load(
            vcf_file=job.lofreq_vcf_path, sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.loaded_at = datetime.now()

    @staticmethod
    def compute_cooccurrence(job: JobEna, queries: Queries, config: Configuration):
        CooccurrenceMatrix().compute(
            sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.cooccurrence_at = datetime.now()
