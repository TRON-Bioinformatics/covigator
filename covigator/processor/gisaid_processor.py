from datetime import datetime

from sqlalchemy import and_
from sqlalchemy.orm import Session

from covigator.misc import backoff_retrier
from covigator.database.model import SampleGisaid, JobStatus, JobGisaid, Sample, DataSource, Log, CovigatorModule
from covigator.database.database import Database, session_scope
from logzero import logger
from dask.distributed import Client
import os
from covigator.processor.downloader import Downloader
from covigator.processor.gisaid_pipeline import Pipeline
from covigator.processor.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class GisaidProcessor:

    def __init__(self, database: Database, dask_client: Client):
        self.start_time = datetime.now()
        self.has_error = False
        self.database = database
        assert self.database is not None, "Empty database"
        self.dask_client = dask_client
        assert self.dask_client is not None, "Empty dask client"

    def process(self):
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
            self.has_error = True
        finally:
            self._write_execution_log(session, count)
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
        with session_scope() as session:
            job = GisaidProcessor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.DOWNLOADED)
            if job is not None:
                try:
                    vcf = GisaidPipeline().run(run_accession=run_accession)
                    logger.info("Processed {}".format(job.run_accession))
                    job.status = JobStatus.PROCESSED
                    job.analysed_at = datetime.now()
                    job.vcf_path = vcf
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.info("Analysis error {} {}".format(run_accession, str(e)))
                    job.status = JobStatus.FAILED_PROCESSING
                    job.failed_at = datetime.now()
                    job.error_message = str(e)
        return run_accession

    @staticmethod
    def load(run_accession: str):
        with session_scope() as session:
            job = GisaidProcessor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.PROCESSED)
            if job is not None:
                try:
                    VcfLoader().load(
                        vcf_file=job.vcf_path,
                        sample=Sample(id=job.run_accession, source=DataSource.GISAID),
                        session=session)
                    logger.info("Loaded {}".format(job.run_accession))
                    job.status = JobStatus.LOADED
                    job.loaded_at = datetime.now()
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.info("Loading error {} {}".format(run_accession, str(e)))
                    job.status = JobStatus.FAILED_LOAD
                    job.failed_at = datetime.now()
                    job.error_message = str(e)
        return run_accession

    @staticmethod
    def find_job_by_accession_and_status(run_accession: str, session: Session, status: JobStatus) -> JobGisaid:
        job = session.query(JobGisaid) \
            .filter(and_(JobGisaid.run_accession == run_accession, JobGisaid.status == status)) \
            .first()
        return job

    @staticmethod
    def find_gisaid_run_by_accession(run_accession: str, session: Session) -> SampleGisaid:
        gisaid_run = session.query(SampleGisaid).filter(SampleGisaid.run_accession == run_accession).first()
        return gisaid_run

    def _write_execution_log(self, session: Session, count):
        end_time = datetime.now()
        session.add(Log(
            start=self.start_time,
            end=end_time,
            source=DataSource.GISAID,
            module=CovigatorModule.PROCESSOR,
            has_error=self.has_error,
            data={
                "processed": count
            }
        ))
        session.commit()

