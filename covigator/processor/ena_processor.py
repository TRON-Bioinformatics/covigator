import json
from datetime import datetime
import pandas as pd
import covigator
from covigator.configuration import Configuration
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorErrorProcessingCoverageResults, CovigatorExcludedSampleBadQualityReads, \
    CovigatorExcludedSampleNarrowCoverage
from covigator.database.model import JobStatus, DataSource, SampleEna
from covigator.database.database import Database
from logzero import logger
from dask.distributed import Client
from covigator.processor.abstract_processor import AbstractProcessor
from covigator.pipeline.ena_pipeline import Pipeline
from covigator.pipeline.vcf_loader import load_vcf


class EnaProcessor(AbstractProcessor):

    def __init__(self, database: Database, dask_client: Client, config: Configuration, wait_time=60):
        logger.info("Initialising ENA processor")
        super().__init__(database, dask_client, DataSource.ENA, config, wait_time=wait_time)

    def _process_run(self, run_accession: str):
        """
        Launches all jobs and returns the futures for the final job only
        """
        # NOTE: here we set the priority of each step to ensure a depth first processing
        future = self.dask_client.submit(EnaProcessor.job, self.config, run_accession, priority=1)
        return future

    @staticmethod
    def job(config: Configuration, run_accession):
        return EnaProcessor.run_job(
            config, run_accession, start_status=JobStatus.QUEUED, end_status=JobStatus.FINISHED,
            error_status=JobStatus.FAILED_PROCESSING, data_source=DataSource.ENA,
            function=EnaProcessor.run_all)

    @staticmethod
    def run_all(sample: SampleEna, queries: Queries, config: Configuration) -> SampleEna:
        sample = EnaProcessor.run_pipeline(sample=sample, queries=queries, config=config)
        sample = EnaProcessor.load(sample=sample, queries=queries, config=config)
        return sample

    @staticmethod
    def run_pipeline(sample: SampleEna, queries: Queries, config: Configuration) -> SampleEna:
        fastq1, fastq2 = sample.get_fastq1_and_fastq2()
        pipeline_result = Pipeline(config=config)\
            .run(run_accession=sample.run_accession, fastq1=fastq1, fastq2=fastq2)
        sample.analysed_at = datetime.now()

        # stores the paths to all files output by pipeline
        sample.lofreq_vcf_path = pipeline_result.lofreq_vcf
        sample.ivar_vcf_path = pipeline_result.ivar_vcf
        sample.gatk_vcf_path = pipeline_result.gatk_vcf
        sample.bcftools_vcf_path = pipeline_result.bcftools_vcf
        sample.lofreq_pangolin_path = pipeline_result.lofreq_pangolin
        sample.ivar_pangolin_path = pipeline_result.ivar_pangolin
        sample.gatk_pangolin_path = pipeline_result.gatk_pangolin
        sample.bcftools_pangolin_path = pipeline_result.bcftools_pangolin
        sample.fastp_path = pipeline_result.fastp_qc
        sample.horizontal_coverage_path = pipeline_result.horizontal_coverage
        sample.vertical_coverage_path = pipeline_result.vertical_coverage

        # stores the covigator version
        sample.covigator_processor_version = covigator.VERSION

        # load FAST JSON into the DB
        sample.qc = json.load(open(pipeline_result.fastp_qc))
        # load horizontal coverage values in the database
        sample = EnaProcessor.load_coverage_results(sample)
        # load pangolin results
        sample = EnaProcessor.load_pangolin(sample=sample, path=sample.lofreq_pangolin_path)

        # NOTE: this is a counterintuititve commit. The VCF loading happening after this may do a legitimate rollback
        # but we don't want to rollback changes in the sample, hence this commit
        queries.session.commit()

        return sample

    @staticmethod
    def load_coverage_results(sample: SampleEna) -> SampleEna:
        try:
            data = pd.read_csv(sample.horizontal_coverage_path,
                               sep="\t",
                               dtype={
                                   'numreads': int,
                                   'covbases': int,
                                   'meandepth': float,
                                   'meanbaseq': float,
                                   'meanmapq': float,
                                   'coverage': float,
                               }
                               )
            data.fillna(0, inplace=True)
            sample.mean_depth = float(data.meandepth.loc[0])
            sample.mean_base_quality = float(data.meanbaseq.loc[0])
            sample.mean_mapping_quality = float(data.meanmapq.loc[0])
            sample.num_reads = int(data.numreads.loc[0])
            sample.covered_bases = int(data.covbases.loc[0])
            sample.coverage = float(data.coverage.loc[0])
        except Exception as e:
            raise CovigatorErrorProcessingCoverageResults(e)

        return sample

    @staticmethod
    def load(sample: SampleEna, queries: Queries, config: Configuration) -> SampleEna:
        if sample.mean_mapping_quality < config.mean_mq_thr or sample.mean_base_quality < config.mean_bq_thr:
            raise CovigatorExcludedSampleBadQualityReads("Mean MQ: {}; mean BCQ: {}".format(
                sample.mean_mapping_quality, sample.mean_base_quality))
        if sample.coverage < config.horizontal_coverage_thr:
            raise CovigatorExcludedSampleNarrowCoverage("Horizontal coverage {} %".format(sample.coverage))
        if not config.skip_vcf_loading:
            load_vcf(
                vcf_file=sample.lofreq_vcf_path, run_accession=sample.run_accession, source=DataSource.ENA, session=queries.session)
            sample.loaded_at = datetime.now()
        return sample
