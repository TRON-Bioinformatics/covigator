from datetime import datetime

from faker import Faker

from covigator.database.queries import Queries
from covigator.misc import backoff_retrier
from covigator.database.model import JobStatus, JobEna, Sample, DataSource
from covigator.database.database import Database, session_scope
from logzero import logger
from dask.distributed import Client
import os
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.processor.cooccurrence_matrix import CooccurrenceMatrix
from covigator.processor.downloader import Downloader
from covigator.processor.ena_pipeline import Pipeline
from covigator.processor.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class DummyProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, batch_size: int):
        logger.info("Initialising dummy processor")
        super().__init__(database, dask_client, DataSource.ENA, batch_size)
        self.database = database
        self.faker = Faker()

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_download = self.dask_client.submit(
            DummyProcessor.run_job, run_accession, JobStatus.QUEUED, JobStatus.DOWNLOADED,
            JobStatus.FAILED_DOWNLOAD, DataSource.ENA, DummyProcessor.download, priority=-1)
        future_process = self.dask_client.submit(
            DummyProcessor.run_job, future_download, JobStatus.DOWNLOADED, JobStatus.PROCESSED,
            JobStatus.FAILED_PROCESSING, DataSource.ENA, DummyProcessor.run_pipeline, priority=1)
        future_delete = self.dask_client.submit(
            DummyProcessor.run_job, future_process, JobStatus.PROCESSED, None,
            None, DataSource.ENA, DummyProcessor.delete, priority=3)
        future_load = self.dask_client.submit(
            DummyProcessor.run_job, future_process, JobStatus.PROCESSED, JobStatus.LOADED,
            JobStatus.FAILED_LOAD, DataSource.ENA, DummyProcessor.load, priority=2)
        future_cooccurrence = self.dask_client.submit(
            DummyProcessor.run_job, future_load, JobStatus.LOADED, JobStatus.COOCCURRENCE,
            JobStatus.FAILED_COOCCURRENCE, DataSource.ENA, DummyProcessor.compute_cooccurrence, priority=2)

        return [future_download, future_process, future_delete, future_load, future_cooccurrence]

    @staticmethod
    def download(job: JobEna, queries: Queries) -> str:
        paths = Faker().filepath(extension="fastq")
        job.fastq_path = paths
        job.downloaded_at = datetime.now()

    @staticmethod
    def run_pipeline(job: JobEna, queries: Queries) -> str:
        job.analysed_at = datetime.now()
        job.vcf_path = Faker().filepath(extension="vcf")

    @staticmethod
    def delete(job: JobEna, queries: Queries):
        # delete FASTQ files from the file system here
        for fastq in job.get_fastq_paths():
            logger.info("Faking deletion of {}".format(fastq))
        job.cleaned_at = datetime.now()

    @staticmethod
    def load(job: JobEna, queries: Queries):
        job.loaded_at = datetime.now()

    @staticmethod
    def compute_cooccurrence(job: JobEna, queries: Queries):
        job.cooccurrence_at = datetime.now()
