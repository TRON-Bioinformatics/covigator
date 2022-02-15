import pandas as pd
from datetime import datetime

from covigator.configuration import Configuration
from covigator.database.model import JobStatus, JobGisaid, Sample, DataSource
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
    def run_all(job: JobGisaid, queries: Queries, config: Configuration):
        GisaidProcessor.run_pipeline(job=job, queries=queries, config=config)
        GisaidProcessor.load(job=job, queries=queries, config=config)

    @staticmethod
    def run_pipeline(job: JobGisaid, queries: Queries, config: Configuration):
        sample = queries.find_sample_by_accession(job.run_accession, source=DataSource.GISAID)
        pipeline_results = GisaidPipeline(config=config).run(sample=sample)
        job.analysed_at = datetime.now()
        job.vcf_path = pipeline_results.vcf_path
        job.pangolin_path = pipeline_results.pangolin_path
        job.fasta_path = pipeline_results.fasta_path

        # load pangolin results
        GisaidProcessor.load_gisaid_pangolin(job)

    @staticmethod
    def load(job: JobGisaid, queries: Queries, config: Configuration):
        VcfLoader().load(
            vcf_file=job.vcf_path, sample=Sample(id=job.run_accession, source=DataSource.GISAID),
            session=queries.session)
        job.loaded_at = datetime.now()

    @staticmethod
    def load_gisaid_pangolin(job: JobGisaid):
        try:
            data = pd.read_csv(job.pangolin_path)
            job.pangolin_lineage = float(data.lineage.loc[0])
            job.pangolin_conflict = float(data.conflict.loc[0])
            job.pangolin_ambiguity_score = float(data.ambiguity_score.loc[0])
            job.pangolin_scorpio_call = float(data.scorpio_call.loc[0])
            job.pangolin_scorpio_support = float(data.scorpio_support.loc[0])
            job.pangolin_scorpio_conflict = float(data.scorpio_conflict.loc[0])
            job.pangolin_version = float(data.version.loc[0])
            job.pangolin_pangolin_version = float(data.pangolin_version.loc[0])
            job.pangolin_pangoLEARN_version = float(data.pangoLEARN_version.loc[0])
            job.pangolin_pango_version = float(data.pango_version.loc[0])
            job.pangolin_status = float(data.status.loc[0])
            job.pangolin_note = float(data.note.loc[0])
        except Exception as e:
            raise CovigatorErrorProcessingPangolinResults(e)
