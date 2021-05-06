from datetime import datetime

from covigator.configuration import Configuration
from covigator.database.model import JobStatus, JobGisaid, Sample, DataSource
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client

from covigator.database.queries import Queries
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
from covigator.pipeline.vcf_loader import VcfLoader

NUMBER_RETRIES_DOWNLOADER = 5


class GisaidProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising GISAID processor")
        super().__init__(database, dask_client, DataSource.GISAID, config=config)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future_process = self.dask_client.submit(
            GisaidProcessor.run_job, self.config, run_accession, JobStatus.QUEUED, JobStatus.PROCESSED,
            JobStatus.FAILED_PROCESSING, DataSource.GISAID, GisaidProcessor.run_pipeline,
            priority=1)
        future_load = self.dask_client.submit(
            GisaidProcessor.run_job, self.config, future_process, JobStatus.PROCESSED, JobStatus.FINISHED,
            JobStatus.FAILED_LOAD, DataSource.GISAID, GisaidProcessor.load,
            priority=2)
        return future_load

    @staticmethod
    def run_pipeline(job: JobGisaid, queries: Queries, config: Configuration):
        sample = queries.find_sample_gisaid_by_accession(job.run_accession)
        vcf = GisaidPipeline(config=config).run(sample=sample)
        job.analysed_at = datetime.now()
        job.vcf_path = vcf
        return job.run_accession

    @staticmethod
    def load(job: JobGisaid, queries: Queries, config: Configuration):
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.GISAID),
            session=queries.session)
        job.loaded_at = datetime.now()
        return job.run_accession
