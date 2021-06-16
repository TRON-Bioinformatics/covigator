import os
import abc
import time
import traceback
from datetime import datetime
from typing import Callable

import typing as typing
from dask.distributed import Client
from sqlalchemy.orm import Session
from logzero import logger
import covigator
import covigator.configuration
from covigator.configuration import Configuration
from covigator.database.database import Database, session_scope
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus, JobEna, JobGisaid
from covigator.database.queries import Queries


class AbstractProcessor:

    def __init__(self, database: Database, dask_client: Client, data_source: DataSource, config: Configuration):
        self.data_source = data_source
        self.config = config
        self.start_time = datetime.now()
        self.has_error = False
        self.error_message = None
        self.database = database
        assert self.database is not None, "Empty database"
        self.dask_client = dask_client
        assert self.dask_client is not None, "Empty dask client"

    def process(self):
        logger.info("Starting processor")
        session = self.database.get_database_session()
        queries = Queries(session)
        count = 0
        try:
            futures = []
            while True:
                job = queries.find_first_pending_job(self.data_source)
                if not job:
                    logger.info("No more jobs to process after sending {} runs to process".format(count))
                    break

                # it has to update the status before doing anything so this processor does not read it again
                job.status = JobStatus.QUEUED
                job.queued_at = datetime.now()
                session.commit()

                # sends the run for processing
                futures.append(self._process_run(run_accession=job.run_accession))
                count += 1
            # waits for all to finish
            self.dask_client.gather(futures=futures)
            logger.info("Processor finished!")
        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.error_message = self._get_traceback_from_exception(e)
            self.has_error = True
        finally:
            self._write_execution_log(session, count, data_source=self.data_source)
            session.close()
            logger.info("Finished processor")

    @staticmethod
    def run_job(config: Configuration, run_accession: str, start_status: JobStatus, end_status: JobStatus,
                error_status: JobStatus, data_source: DataSource,
                function: Callable[[typing.Union[JobEna, JobGisaid], Queries, Configuration], None]) -> str or None:
        """
        Runs a function on a job, if anything goes wrong or does not fit in the DB it returns None in order to
        stop the execution of subsequent jobs.
        """
        covigator.configuration.initialise_logs(config.logfile_processor, sample_id=run_accession)
        if run_accession is not None:
            try:
                with session_scope(config=config) as session:
                    queries = Queries(session)
                    job = queries.find_job_by_accession_and_status(
                        run_accession=run_accession, status=start_status, data_source=data_source)
                    if job is not None:
                        function(job, queries, config)
                        if end_status is not None:
                            job.status = end_status
                    else:
                        logger.warning("Expected ENA job {} in status {}".format(run_accession, start_status))
                        run_accession = None
            except Exception as e:
                # captures any possible exception happening, but logs it in the DB
                if error_status is not None:
                    AbstractProcessor._log_error_in_job(
                        config=config, run_accession=run_accession, exception=e, status=error_status, data_source=data_source)
                    run_accession = None
                else:
                    logger.warning("Error processing a job that does not stop the workflow!")
        return run_accession

    @abc.abstractmethod
    def _process_run(self, run_accession: str):
        pass

    def _write_execution_log(self, session: Session, count, data_source: DataSource):
        end_time = datetime.now()
        session.add(Log(
            start=self.start_time,
            end=end_time,
            source=data_source,
            module=CovigatorModule.PROCESSOR,
            has_error=self.has_error,
            processed=count,
            data={
                "processed": count
            }
        ))
        session.commit()

    @staticmethod
    def _get_traceback_from_exception(e):
        return "".join(traceback.format_exception(
            etype=type(e), value=e, tb=e.__traceback__))

    @staticmethod
    def _log_error_in_job(config: Configuration, run_accession: str, exception: Exception, status: JobStatus, data_source: DataSource):
        with session_scope(config=config) as session:
            logger.info("Error on job {} on state {}: {}".format(run_accession, status, str(exception)))
            job = Queries(session).find_job_by_accession(run_accession=run_accession, data_source=data_source)
            job.status = status
            job.failed_at = datetime.now()
            job.error_message = AbstractProcessor._get_traceback_from_exception(exception)