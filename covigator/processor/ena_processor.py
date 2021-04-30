import uuid
from datetime import datetime

from covigator.configuration import Configuration
from covigator.database.queries import Queries
from covigator.misc import backoff_retrier
from covigator.database.model import JobStatus, JobEna, Sample, DataSource
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client
import os
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.cooccurrence_matrix import CooccurrenceMatrix
from covigator.pipeline.downloader import Downloader
from covigator.pipeline.ena_pipeline import Pipeline
from covigator.pipeline.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class EnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_download = self.dask_client.submit(EnaProcessor.download_job, self.config, run_accession, priority=-1)
        future_process = self.dask_client.submit(EnaProcessor.pipeline_job, self.config, future_download, priority=1)
        future_delete = self.dask_client.submit(EnaProcessor.cleanup_job, self.config, future_process, priority=3)
        future_load = self.dask_client.submit(EnaProcessor.load_job, self.config, future_process, priority=2)
        future_cooccurrence = self.dask_client.submit(
            EnaProcessor.cooccurrence_job, self.config, future_load, priority=2)

        return [future_download, future_process, future_delete, future_load, future_cooccurrence]

    @staticmethod
    def cooccurrence_job(config: Configuration, run_accession):
        EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.LOADED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_COOCCURRENCE, data_source=DataSource.ENA,
            function=EnaProcessor.compute_cooccurrence)

    @staticmethod
    def load_job(config: Configuration, run_accession):
        EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.PROCESSED, end_status=JobStatus.LOADED,
            error_status=JobStatus.FAILED_LOAD, data_source=DataSource.ENA, function=EnaProcessor.load)

    @staticmethod
    def download_job(config: Configuration, run_accession):
        EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.DOWNLOADED,
            error_status=JobStatus.FAILED_DOWNLOAD, data_source=DataSource.ENA, function=EnaProcessor.download)

    @staticmethod
    def pipeline_job(config: Configuration, run_accession):
        EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.DOWNLOADED, end_status=JobStatus.PROCESSED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA, function=EnaProcessor.run_pipeline)

    @staticmethod
    def cleanup_job(config: Configuration, run_accession):
        EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.PROCESSED, end_status=None, error_status=None,
            data_source=DataSource.ENA, function=EnaProcessor.delete)

    @staticmethod
    def download(job: JobEna, queries: Queries, config: Configuration):
        # ensures that the download is done with retries, even after MD5 check sum failure
        download_with_retries = backoff_retrier.wrapper(Downloader(config=config).download, NUMBER_RETRIES_DOWNLOADER)
        sample_ena = queries.find_sample_ena_by_accession(job.run_accession)
        paths = download_with_retries(sample_ena=sample_ena)
        job.fastq_path = paths
        job.downloaded_at = datetime.now()

    @staticmethod
    def run_pipeline(job: JobEna, queries: Queries, config: Configuration):
        fastq1, fastq2 = job.get_fastq1_and_fastq2()
        vcf = Pipeline(config=config).run(fastq1=fastq1, fastq2=fastq2)
        job.analysed_at = datetime.now()
        job.vcf_path = vcf

    @staticmethod
    def delete(job: JobEna, queries: Queries, config: Configuration):
        # delete FASTQ files from the file system here
        for fastq in job.get_fastq_paths():
            os.remove(fastq)
        job.cleaned_at = datetime.now()

    @staticmethod
    def load(job: JobEna, queries: Queries, config: Configuration):
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.loaded_at = datetime.now()

    @staticmethod
    def compute_cooccurrence(job: JobEna, queries: Queries, config: Configuration):
        CooccurrenceMatrix().compute(
            sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.cooccurrence_at = datetime.now()
