from datetime import datetime

from sqlalchemy import and_
from sqlalchemy.orm import Session

from covigator.misc import backoff_retrier
from covigator.model import Database, EnaRun, JobStatus, Job
from logzero import logger
from dask.distributed import Client
import os
from covigator.processor.downloader import Downloader
from covigator.processor.pipeline import Pipeline

NUMBER_RETRIES_DOWNLOADER = 5


class Processor:

    def __init__(self, database: Database, dask_client: Client):
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
                job = session.query(Job)\
                    .filter(Job.status == JobStatus.PENDING)\
                    .order_by(Job.created_at.desc())\
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
        finally:
            session.close()

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_download = self.dask_client.submit(
            Processor.download, run_accession, priority=-1)
        future_process = self.dask_client.submit(
            Processor.run_pipeline, future_download, priority=1)
        future_delete = self.dask_client.submit(
            Processor.delete, future_process, priority=3)
        future_load = self.dask_client.submit(
            Processor.load, future_process, priority=2)
        return [future_download, future_process, future_delete, future_load]

    @staticmethod
    def download(run_accession: str) -> str:
        session = Database().get_database_session()
        try:
            job = Processor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.QUEUED)
            # NOTE: eventually we may want to generalize this to support some other jobs than ENA runs
            ena_run = Processor.find_ena_run_by_accession(run_accession, session)
            if job is not None:
                try:
                    # ensures that the download is done with retries, even after MD5 check sum failure
                    download_with_retries = backoff_retrier.wrapper(Downloader().download, NUMBER_RETRIES_DOWNLOADER)
                    paths = download_with_retries(ena_run=ena_run)
                    job.fastq_path = paths
                    job.status = JobStatus.DOWNLOADED
                    job.downloaded_at = datetime.now()
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.info("Download error {} {}".format(run_accession, str(e)))
                    job.status = JobStatus.FAILED_DOWNLOAD
                    job.failed_at = datetime.now()
                    job.error_message = str(e)
            session.commit()
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()
        return run_accession

    @staticmethod
    def run_pipeline(run_accession: str) -> str:
        session = Database().get_database_session()
        try:
            job = Processor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.DOWNLOADED)
            if job is not None:
                try:
                    # TODO: get the path to the VCF coming out of the pipeline
                    Pipeline(fastqs=job.get_fastq_paths()).run()
                    logger.info("Processed {}".format(job.run_accession))
                    job.status = JobStatus.PROCESSED
                    job.analysed_at = datetime.now()
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.info("Analysis error {} {}".format(run_accession, str(e)))
                    job.status = JobStatus.FAILED_PROCESSING
                    job.failed_at = datetime.now()
                    job.error_message = str(e)
            session.commit()
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()
        return run_accession

    @staticmethod
    def delete(run_accession: str):
        session = Database().get_database_session()
        try:
            job = Processor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.PROCESSED)
            if job is not None:
                try:
                    # delete FASTQ files from the file system here
                    for fastq in job.get_fastq_paths():
                        os.remove(fastq)
                    logger.info("Deleted {}".format(job.run_accession))
                    job.cleaned_at = datetime.now()
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.error("File deletion error {} {}".format(run_accession, str(e)))
                    # we don't do anything as cleaning up is not really a must in the pipeline
            session.commit()
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()
        return run_accession

    @staticmethod
    def load(run_accession: str):
        session = Database().get_database_session()
        try:
            job = Processor.find_job_by_accession_and_status(
                run_accession=run_accession, session=session, status=JobStatus.PROCESSED)
            if job is not None:
                try:
                    # TODO: load VCFs in DB
                    logger.info("Loaded {}".format(job.run_accession))
                    job.status = JobStatus.LOADED
                    job.loaded_at = datetime.now()
                except Exception as e:  # TODO: do we want a less wide exception capture?
                    logger.info("Loading error {} {}".format(run_accession, str(e)))
                    job.status = JobStatus.FAILED_LOAD
                    job.failed_at = datetime.now()
                    job.error_message = str(e)
            session.commit()
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()
        return run_accession

    @staticmethod
    def find_job_by_accession_and_status(run_accession: str, session: Session, status: JobStatus) -> Job:
        job = session.query(Job) \
            .filter(and_(Job.run_accession == run_accession, Job.status == status)) \
            .first()
        return job

    @staticmethod
    def find_ena_run_by_accession(run_accession: str, session: Session) -> EnaRun:
        ena_run = session.query(EnaRun).filter(EnaRun.run_accession == run_accession).first()
        return ena_run
