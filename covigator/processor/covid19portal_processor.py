from datetime import datetime

import covigator
from covigator.configuration import Configuration
from covigator.database.model import JobStatus, DataSource, SampleCovid19Portal
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client

from covigator.database.queries import Queries
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.covid19_portal_pipeline import Covid19PortalPipeline
from covigator.pipeline.vcf_loader import VcfLoader


class Covid19PortalProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration, wait_time=60):
        logger.info("Initialising Covid19 Portal processor")
        super().__init__(database, dask_client, DataSource.COVID19_PORTAL, config=config, wait_time=wait_time)

    def _process_run(self, run_accession: str):
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(Covid19PortalProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return Covid19PortalProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.COVID19_PORTAL,
            function=Covid19PortalProcessor.run_all)

    @staticmethod
    def run_all(sample: SampleCovid19Portal, queries: Queries, config: Configuration):
        sample = Covid19PortalProcessor.run_pipeline(sample=sample, queries=queries, config=config)
        sample = Covid19PortalProcessor.load(sample=sample, queries=queries, config=config)
        return sample

    @staticmethod
    def run_pipeline(sample: SampleCovid19Portal, queries: Queries, config: Configuration) -> SampleCovid19Portal:

        pipeline_results = Covid19PortalPipeline(config=config).run(sample=sample)
        sample.analysed_at = datetime.now()
        sample.sample_folder = sample.get_sample_folder(config.storage_folder)
        sample.vcf_path = pipeline_results.vcf_path
        sample.pangolin_path = pipeline_results.pangolin_path
        sample.fasta_path = pipeline_results.fasta_path
        sample.covigator_processor_version = covigator.VERSION  # stores the covigator version

        # NOTE: this is a counterintuititve commit. The VCF loading happening after this may do a legitimate rollback
        # but we don't want to rollback changes in the sample, hence this commit
        queries.session.commit()

        return sample

    @staticmethod
    def load(sample: SampleCovid19Portal, queries: Queries, config: Configuration) -> SampleCovid19Portal:
        if not config.skip_vcf_loading:
            VcfLoader().load(
                vcf_file=sample.vcf_path, run_accession=sample.run_accession, source=DataSource.COVID19_PORTAL,
                session=queries.session)
            sample.loaded_at = datetime.now()
        return sample
