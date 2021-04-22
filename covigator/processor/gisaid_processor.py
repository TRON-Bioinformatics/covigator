from datetime import datetime
from covigator.database.model import JobStatus, JobGisaid, Sample, DataSource
from covigator.database.database import Database, session_scope
from logzero import logger
from dask.distributed import Client

from covigator.database.queries import Queries
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.processor.gisaid_pipeline import GisaidPipeline
from covigator.processor.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class GisaidProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client):
        logger.info("Initialising GISAID processor")
        super().__init__(database, dask_client, DataSource.GISAID)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_process = self.dask_client.submit(
            GisaidProcessor.run_pipeline, run_accession, priority=1)
        future_load = self.dask_client.submit(
            GisaidProcessor.load, future_process, priority=2)
        return [future_process, future_load]

    @staticmethod
    def run_pipeline(run_accession: str) -> str:
        try:
            with session_scope() as session:
                queries = Queries(session)
                job = queries.find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.QUEUED, data_source=DataSource.GISAID)
                if job is not None:
                    vcf = GisaidPipeline().run(run_accession=run_accession)
                    logger.info("Processed {}".format(job.run_accession))
                    job.status = JobStatus.PROCESSED
                    job.analysed_at = datetime.now()
                    job.vcf_path = vcf
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            GisaidProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_PROCESSING, data_source=DataSource.GISAID)
        return run_accession

    @staticmethod
    def load(run_accession: str):
        try:
            with session_scope() as session:
                queries = Queries(session)
                job = queries.find_job_by_accession_and_status(
                    run_accession=run_accession, status=JobStatus.PROCESSED, data_source=DataSource.GISAID)
                if job is not None:
                    VcfLoader().load(
                        vcf_file=job.vcf_path,
                        sample=Sample(id=job.run_accession, source=DataSource.GISAID),
                        session=session)
                    logger.info("Loaded {}".format(job.run_accession))
                    job.status = JobStatus.LOADED
                    job.loaded_at = datetime.now()
        except Exception as e:
            # captures any possible exception happening, but logs it in the DB
            GisaidProcessor._log_error_in_job(
                run_accession=run_accession, exception=e, status=JobStatus.FAILED_LOAD, data_source=DataSource.GISAID)
        return run_accession
