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
        super().__init__(database, dask_client)

    def process(self):
        logger.info("Starting GISAID processor")
        session = self.database.get_database_session()
        count = 0
        try:
            futures = []
            while True:
                job = session.query(JobGisaid)\
                    .filter(JobGisaid.status == JobStatus.PENDING)\
                    .order_by(JobGisaid.created_at.desc())\
                    .first()
                if not job:
                    logger.info("No more jobs to process after sending {} runs to process".format(count))
                    break

                # it has to update the status before doing anything so this processor does not read it again
                job.status = JobStatus.QUEUED
                job.queued_at = datetime.now()
                session.commit()

                # sends the run for processing
                futures.extend(self._process_run(run_accession=job.run_accession))
                count += 1
            self.dask_client.gather(futures=futures)

        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.error_message = self._get_traceback_from_exception(e)
            self.has_error = True
        finally:
            self._write_execution_log(session, count, data_source=DataSource.GISAID)
            session.close()

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
