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


class GisaidProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration):
        logger.info("Initialising GISAID processor")
        super().__init__(database, dask_client, DataSource.GISAID, config=config)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(GisaidProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return GisaidProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.GISAID,
            function=GisaidProcessor.run_all)

    @staticmethod
    def run_all(job: JobGisaid, queries: Queries, config: Configuration):
        GisaidProcessor.run_pipeline(job=job, queries=queries, config=config)
        GisaidProcessor.load(job=job, queries=queries, config=config)

    @staticmethod
    def run_pipeline(job: JobGisaid, queries: Queries, config: Configuration):
        sample = queries.find_sample_by_accession(job.run_accession, source=DataSource.GISAID)
        vcf = GisaidPipeline(config=config).run(sample=sample)
        job.analysed_at = datetime.now()
        job.vcf_path = vcf

    @staticmethod
    def load(job: JobGisaid, queries: Queries, config: Configuration):
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.GISAID),
            session=queries.session,
            max_snvs=config.max_snvs, max_deletions=config.max_deletions, max_insertions=config.max_insertions)
        job.loaded_at = datetime.now()
