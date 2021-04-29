from datetime import datetime
from covigator.database.queries import Queries
from covigator.misc import backoff_retrier
from covigator.database.model import JobStatus, JobEna, Sample, DataSource
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client
import os
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.processor.cooccurrence_matrix import CooccurrenceMatrix
from covigator.processor.downloader import Downloader
from covigator.processor.ena_pipeline import Pipeline
from covigator.processor.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class EnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_download = self.dask_client.submit(
            EnaProcessor.run_job, run_accession, JobStatus.QUEUED, JobStatus.DOWNLOADED,
            JobStatus.FAILED_DOWNLOAD, DataSource.ENA, EnaProcessor.download, priority=-1)
        future_process = self.dask_client.submit(
            EnaProcessor.run_job, future_download, JobStatus.DOWNLOADED, JobStatus.PROCESSED,
            JobStatus.FAILED_PROCESSING, DataSource.ENA, EnaProcessor.run_pipeline, priority=1)
        future_delete = self.dask_client.submit(
            EnaProcessor.run_job, future_process, JobStatus.PROCESSED, None,
            None, DataSource.ENA, EnaProcessor.delete, priority=3)
        future_load = self.dask_client.submit(
            EnaProcessor.run_job, future_process, JobStatus.PROCESSED, JobStatus.LOADED,
            JobStatus.FAILED_LOAD, DataSource.ENA, EnaProcessor.load, priority=2)
        future_cooccurrence = self.dask_client.submit(
            EnaProcessor.run_job, future_load, JobStatus.LOADED, JobStatus.FINISHED,
            JobStatus.FAILED_COOCCURRENCE, DataSource.ENA, EnaProcessor.compute_cooccurrence, priority=2)

        return [future_download, future_process, future_delete, future_load, future_cooccurrence]

    @staticmethod
    def download(job: JobEna, queries: Queries):
        # ensures that the download is done with retries, even after MD5 check sum failure
        download_with_retries = backoff_retrier.wrapper(Downloader().download, NUMBER_RETRIES_DOWNLOADER)
        sample_ena = queries.find_sample_ena_by_accession(job.run_accession)
        paths = download_with_retries(sample_ena=sample_ena)
        job.fastq_path = paths
        job.downloaded_at = datetime.now()

    @staticmethod
    def run_pipeline(job: JobEna, queries: Queries):
        fastq1, fastq2 = job.get_fastq1_and_fastq2()
        vcf = Pipeline().run(fastq1=fastq1, fastq2=fastq2)
        job.analysed_at = datetime.now()
        job.vcf_path = vcf

    @staticmethod
    def delete(job: JobEna, queries: Queries):
        # delete FASTQ files from the file system here
        for fastq in job.get_fastq_paths():
            os.remove(fastq)
        job.cleaned_at = datetime.now()

    @staticmethod
    def load(job: JobEna, queries: Queries):
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.loaded_at = datetime.now()

    @staticmethod
    def compute_cooccurrence(job: JobEna, queries: Queries):
        CooccurrenceMatrix().compute(
            sample=Sample(id=job.run_accession, source=DataSource.ENA), session=queries.session)
        job.cooccurrence_at = datetime.now()
