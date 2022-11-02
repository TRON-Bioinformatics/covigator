import abc
import time
import traceback
import pandas as pd
from contextlib import suppress
from datetime import datetime, date
from typing import Callable
import typing as typing
from dask.distributed import Client
from distributed import fire_and_forget, wait
from logzero import logger
from sqlalchemy.exc import SQLAlchemyError

from covigator.configuration import Configuration
from covigator.database.database import Database, session_scope
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus, \
    SampleEna, LastUpdate
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorExcludedSampleException, CovigatorErrorProcessingPangolinResults
from covigator.precomputations.loader import PrecomputationsLoader


class AbstractProcessor:

    def __init__(self, database: Database, dask_client: Client, data_source: DataSource, config: Configuration,
                 wait_time=60):
        self.data_source = data_source
        self.config = config
        self.start_time = datetime.now()
        self.has_error = False
        self.error_message = None
        self.database = database
        assert self.database is not None, "Empty database"
        self.dask_client = dask_client
        assert self.dask_client is not None, "Empty dask client"
        self.session = self.database.get_database_session()
        self.queries = Queries(self.session)
        self.wait_time = wait_time

    def process(self):
        logger.info("Starting processor")
        count = 0
        try:
            futures = []
            while True:
                # queries 100 jobs every time to make sending to queue faster
                jobs = self.queries.find_first_pending_jobs(self.data_source, n=1000, status=[JobStatus.DOWNLOADED])
                if jobs is None or len(jobs) == 0:
                    logger.info("No more jobs to process after sending {} runs to process".format(count))
                    break
                for job in jobs:
                    # it has to update the status before doing anything so a processor does not read it again
                    try:
                        job.status = JobStatus.QUEUED
                        job.queued_at = datetime.now()
                        self.session.commit()
                    except SQLAlchemyError:
                        # this job has been taken by another processor
                        continue

                    # sends the run for processing
                    future = self._process_run(run_accession=job.run_accession)
                    futures.extend(future)
                    count += 1
                    if count % 1000 == 0:
                        logger.info("Sent {} jobs for processing...".format(count))

                    # waits for a batch to finish
                    if len(futures) >= self.config.batch_size:
                        # waits for a batch to finish before sending more
                        logger.info("Waiting for a batch to be processed...")
                        wait(fs=futures)
                        futures = []
                        logger.info("Batch finished!")

            # waits for the last batch to finish
            if len(futures) > 0:
                logger.info("Waiting for the last batch to be processed...")
                wait(fs=futures)
            logger.info("Processor finished!")

            # precomputes data right after processor
            PrecomputationsLoader(session=self.session).load()

            # updates the last update entry
            self._register_last_update()

        except Exception as e:
            logger.exception(e)
            self.session.rollback()
            self.error_message = self._get_traceback_from_exception(e)
            self.has_error = True
        finally:
            logger.info("Logging execution stats...")
            self._write_execution_log(count, data_source=self.data_source)
            logger.info("Waits {} secs to let the cluster tidy up things...".format(self.wait_time))
            time.sleep(self.wait_time)
            logger.info("Shutting down cluster and database session...")
            with suppress(Exception):
                self.dask_client.close(60)
                self.session.close()
            logger.info("Cluster and database sessions closed")

    def _register_last_update(self):
        last_update = LastUpdate(source=self.data_source, update_time=date.today())
        self.session.add(last_update)
        self.session.commit()

    @staticmethod
    def run_job(config: Configuration, run_accession: str, start_status: JobStatus, end_status: JobStatus,
                error_status: JobStatus, data_source: DataSource,
                function: Callable[[typing.Union[SampleEna], Queries, Configuration],
                                   typing.Union[SampleEna]]) -> str or None:
        """
        Runs a function on a job, if anything goes wrong or does not fit in the DB it returns None in order to
        stop the execution of subsequent jobs.
        """
        if run_accession is not None:
            try:
                with session_scope(config=config) as session:
                    queries = Queries(session)
                    sample = queries.find_job_by_accession_and_status(
                        run_accession=run_accession, status=start_status, data_source=data_source)
                    if sample is not None:
                        sample = function(sample, queries, config)
                        if end_status is not None:
                            sample.status = end_status
                    else:
                        run_accession = None
            except CovigatorExcludedSampleException as e:
                # captures exclusion cases
                AbstractProcessor._log_error_in_job(
                    config=config, run_accession=run_accession, exception=e, status=JobStatus.EXCLUDED,
                    data_source=data_source)
                run_accession = None
            except Exception as e:
                # captures any possible exception happening, but logs it in the DB
                if error_status is not None:
                    AbstractProcessor._log_error_in_job(
                        config=config, run_accession=run_accession, exception=e, status=error_status,
                        data_source=data_source)
                    run_accession = None
        return run_accession

    @abc.abstractmethod
    def _process_run(self, run_accession: str):
        pass

    def _write_execution_log(self, count, data_source: DataSource):
        end_time = datetime.now()
        self.session.add(Log(
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
        self.session.commit()

    @staticmethod
    def _get_traceback_from_exception(e):
        return "".join(traceback.format_exception(
            etype=type(e), value=e, tb=e.__traceback__))

    @staticmethod
    def _log_error_in_job(config: Configuration, run_accession: str, exception: Exception, status: JobStatus, data_source: DataSource):
        with session_scope(config=config) as session:
            sample = Queries(session).find_job_by_accession(run_accession=run_accession, data_source=data_source)
            sample.status = status
            sample.failed_at = datetime.now()
            sample.error_message = AbstractProcessor._get_traceback_from_exception(exception)

    @staticmethod
    def load_pangolin(sample: typing.Union[SampleEna], path: str) -> typing.Union[SampleEna]:
        try:
            data = pd.read_csv(path,
                               na_values=None,
                               dtype={
                                   'lineage': str,
                                   'conflict': float,
                                   'ambiguity_score': float,
                                   'scorpio_call': str,
                                   'scorpio_support': float,
                                   'scorpio_conflict': float,
                                   'version': str,
                                   'pangolin_version': str,
                                   'scorpio_version': str,
                                   'constellation_version': str,
                                   'qc_status': str,
                                   'qc_notes': str,
                                   'note': str
                               })

            # fill NA values on a per column basis...
            data.lineage.fillna(value="", inplace=True)
            data.lineage = data.lineage.transform(lambda x: "" if x == 'None' else x)   # replace "None" values
            data.scorpio_call.fillna(value="", inplace=True)
            data.version.fillna(value="", inplace=True)
            data.pangolin_version.fillna(value="", inplace=True)
            data.scorpio_version.fillna(value="", inplace=True)
            data.constellation_version.fillna(value="", inplace=True)
            data.note.fillna(value="", inplace=True)
            data.conflict.fillna(value=0.0, inplace=True)
            data.ambiguity_score.fillna(value=0.0, inplace=True)
            data.scorpio_support.fillna(value=0.0, inplace=True)
            data.scorpio_conflict.fillna(value=0.0, inplace=True)

            sample.pangolin_lineage = data.lineage.loc[0]
            sample.pangolin_conflict = data.conflict.loc[0]
            sample.pangolin_ambiguity_score = data.ambiguity_score.loc[0]
            sample.pangolin_scorpio_call = data.scorpio_call.loc[0]
            sample.pangolin_scorpio_support = data.scorpio_support.loc[0]
            sample.pangolin_scorpio_conflict = data.scorpio_conflict.loc[0]
            sample.pangolin_version = data.version.loc[0]
            sample.pangolin_pangolin_version = data.pangolin_version.loc[0]
            sample.pangolin_scorpio_version = data.scorpio_version.loc[0]
            sample.pangolin_constellation_version = data.constellation_version.loc[0]
            sample.pangolin_qc_status = data.qc_status.loc[0]
            sample.pangolin_qc_notes = data.qc_notes.loc[0]
            sample.pangolin_note = data.note.loc[0]
        except Exception as e:
            raise CovigatorErrorProcessingPangolinResults(e)

        return sample