import abc
import traceback
from datetime import datetime
from dask.distributed import Client
from sqlalchemy.orm import Session
from logzero import logger
from covigator.database.database import Database, session_scope
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus
from covigator.database.queries import Queries


class AbstractProcessor:

    def __init__(self, database: Database, dask_client: Client, data_source: DataSource, batch_size: int):
        self.data_source = data_source
        self.batch_size = batch_size
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
                futures.extend(self._process_run(run_accession=job.run_accession))
                count += 1
                if count % self.batch_size == 0:
                    # process a batch
                    self.dask_client.gather(futures=futures)
                    futures = []
            # process the remaining batch
            self.dask_client.gather(futures=futures)

        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.error_message = self._get_traceback_from_exception(e)
            self.has_error = True
        finally:
            self._write_execution_log(session, count, data_source=self.data_source)
            session.close()
            logger.info("Finished processor")

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
    def _log_error_in_job(run_accession: str, exception: Exception, status: JobStatus, data_source: DataSource):
        with session_scope() as session:
            logger.info("Error on job {} on state {}: {}".format(run_accession, status, str(exception)))
            job = Queries(session).find_job_by_accession(run_accession=run_accession, data_source=data_source)
            job.status = status
            job.failed_at = datetime.now()
            job.error_message = AbstractProcessor._get_traceback_from_exception(exception)