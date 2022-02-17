import pandas as pd
from datetime import datetime

from covigator.configuration import Configuration
from covigator.database.model import JobStatus, DataSource, SampleGisaid
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client

from covigator.database.queries import Queries
from covigator.exceptions import CovigatorErrorProcessingPangolinResults
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
from covigator.pipeline.vcf_loader import VcfLoader


class GisaidProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration, wait_time=60):
        logger.info("Initialising GISAID processor")
        super().__init__(database, dask_client, DataSource.GISAID, config=config, wait_time=wait_time)

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
    def run_all(sample: SampleGisaid, queries: Queries, config: Configuration):
        GisaidProcessor.run_pipeline(sample=sample, queries=queries, config=config)
        GisaidProcessor.load(sample=sample, queries=queries, config=config)

    @staticmethod
    def run_pipeline(sample: SampleGisaid, queries: Queries, config: Configuration):
        pipeline_results = GisaidPipeline(config=config).run(sample=sample)
        sample.analysed_at = datetime.now()
        sample.vcf_path = pipeline_results.vcf_path
        sample.pangolin_path = pipeline_results.pangolin_path
        sample.fasta_path = pipeline_results.fasta_path

        # load pangolin results
        GisaidProcessor.load_pangolin(sample=sample, path=sample.pangolin_path)

    @staticmethod
    def load(sample: SampleGisaid, queries: Queries, config: Configuration):
        VcfLoader().load(
            vcf_file=sample.vcf_path, run_accession=sample.run_accession, source=DataSource.GISAID, session=queries.session)
        sample.loaded_at = datetime.now()
