import os
import pkg_resources
from dask.distributed import Client
from logzero import logger

import covigator
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.database.model import DataSource, JobStatus, JobEna, JobGisaid
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorExcludedSampleTooManyMutations
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.processor.ena_processor import EnaProcessor
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID


class FakeConfiguration(Configuration):

    def __init__(self):
        # use a folder with write permissions for testing
        os.environ[self.ENV_COVIGATOR_STORAGE_FOLDER] = "./data/covigator"
        os.environ[self.ENV_COVIGATOR_REF_FASTA] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/MN908947.3.fa")

        # this makes sure that we do not wipe a relevant database by mistake
        os.environ[self.ENV_COVIGATOR_TABLE_VERSION] = "_test"

        super().__init__()


class FakeEnaAccessor(EnaAccessor):

    def __init__(self, results, database=None):
        # uses an in memory database or the one provided
        super().__init__(tax_id=SARS_COV_2_TAXID, host_tax_id=HOMO_SAPIENS_TAXID,
                         database=database if database else Database(test=True, config=Configuration()))
        self.results = results

    def _get_ena_runs_page(self):
        return self.results


class FakeEnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(FakeEnaProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return FakeEnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA,
            function=FakeEnaProcessor.run_all
    )

    @staticmethod
    def run_all(job: JobEna, queries: Queries, config: Configuration):
        logger.info("Job processed!")


class FakeEnaProcessorExcludingSamples(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(FakeEnaProcessorExcludingSamples.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return FakeEnaProcessorExcludingSamples.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA,
            function=FakeEnaProcessorExcludingSamples.run_all
        )

    @staticmethod
    def run_all(job: JobEna, queries: Queries, config: Configuration):
        raise CovigatorExcludedSampleTooManyMutations("Exclude em'all")


class FakeEnaProcessorFailing(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(FakeEnaProcessorFailing.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return FakeEnaProcessorFailing.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA,
            function=FakeEnaProcessorFailing.run_all
        )

    @staticmethod
    def run_all(job: JobEna, queries: Queries, config: Configuration):
        raise ValueError("Fail em'all")


class FakeGisaidProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising GISAID processor")
        super().__init__(database, dask_client, DataSource.GISAID, config)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(FakeGisaidProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return FakeGisaidProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.GISAID,
            function=FakeGisaidProcessor.run_all
    )

    @staticmethod
    def run_all(job: JobGisaid, queries: Queries, config: Configuration):
        logger.info("Job processed!")
