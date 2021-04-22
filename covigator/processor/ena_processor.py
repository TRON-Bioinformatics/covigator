from datetime import datetime
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


class EnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_download = self.dask_client.submit(
            EnaProcessor.download, run_accession, priority=-1)
        future_process = self.dask_client.submit(
            EnaProcessor.run_pipeline, future_download, priority=1)
        future_delete = self.dask_client.submit(
            EnaProcessor.delete, future_process, priority=3)
        future_load = self.dask_client.submit(
            EnaProcessor.load, future_process, priority=2)
        future_cooccurrence = self.dask_client.submit(
            EnaProcessor.compute_cooccurrence, future_load, priority=2)
        return [future_download, future_process, future_delete, future_load, future_cooccurrence]

    @staticmethod
    def download(run_accession: str) -> str:
        try:
            with session_scope() as session:
                queries = Queries(session)
                job = queries.find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.QUEUED, data_source=DataSource.GISAID)
                # NOTE: eventually we may want to generalize this to support some other jobs than ENA runs
                ena_run = queries.find_ena_run_by_accession(run_accession)
                if job is not None:
                    # ensures that the download is done with retries, even after MD5 check sum failure
                    download_with_retries = backoff_retrier.wrapper(Downloader().download, NUMBER_RETRIES_DOWNLOADER)
                    paths = download_with_retries(ena_run=ena_run)
                    job.fastq_path = paths
                    job.status = JobStatus.DOWNLOADED
                    job.downloaded_at = datetime.now()
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            EnaProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_DOWNLOAD,
                data_source=DataSource.ENA)
        return run_accession

    @staticmethod
    def run_pipeline(run_accession: str) -> str:
        try:
            with session_scope() as session:
                job = Queries(session).find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.DOWNLOADED, data_source=DataSource.ENA)
                if job is not None:
                    fastq1, fastq2 = job.get_fastq1_and_fastq2()
                    vcf = Pipeline().run(fastq1=fastq1, fastq2=fastq2)
                    logger.info("Processed {}".format(job.run_accession))
                    job.status = JobStatus.PROCESSED
                    job.analysed_at = datetime.now()
                    job.vcf_path = vcf
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            EnaProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_PROCESSING,
                data_source=DataSource.ENA)
        return run_accession

    @staticmethod
    def delete(run_accession: str):
        try:
            with session_scope() as session:
                job = Queries(session).find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.PROCESSED, data_source=DataSource.ENA)
                if job is not None:
                    # delete FASTQ files from the file system here
                    for fastq in job.get_fastq_paths():
                        os.remove(fastq)
                    logger.info("Deleted {}".format(job.run_accession))
                    job.cleaned_at = datetime.now()
        except Exception as e:
            # we don't do track in the db if cleanup fails
            logger.warning("Clean up for job {} failed: {}".format(run_accession, str(e)))
        return run_accession

    @staticmethod
    def load(run_accession: str):
        try:
            with session_scope() as session:
                job = Queries(session).find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.PROCESSED, data_source=DataSource.ENA)
                if job is not None:
                    VcfLoader().load(
                        vcf_file=job.vcf_path,
                        sample=Sample(id=job.run_accession, source=DataSource.ENA),
                        session=session)
                    logger.info("Loaded {}".format(job.run_accession))
                    job.status = JobStatus.LOADED
                    job.loaded_at = datetime.now()
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            EnaProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_LOAD, data_source=DataSource.ENA)
        return run_accession

    @staticmethod
    def compute_cooccurrence(run_accession: str):
        try:
            with session_scope() as session:
                job = Queries(session).find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.LOADED, data_source=DataSource.ENA)
                if job is not None:
                    CooccurrenceMatrix().compute(
                        sample=Sample(id=job.run_accession, source=DataSource.ENA),
                        session=session)
                    logger.info("Cooccurrence matrix computed {}".format(job.run_accession))
                    job.status = JobStatus.COOCCURRENCE
                    job.cooccurrence_at = datetime.now()
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            EnaProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_COOCCURRENCE,
                data_source=DataSource.ENA)
        return run_accession
