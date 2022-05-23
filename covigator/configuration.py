import logging
import os
import logzero
from logzero import logger

import covigator
from covigator.exceptions import CovigatorDashBoardInitialisationError

DEFAULT_SUBCLONAL_AF_THR = 0.8
DEFAULT_LOW_COVERAGE_AF_THR = 0.05
DEFAULT_HORIZONTAL_COVERAGE_THR = 20
DEFAULT_MEAN_BQ_THR = 10
DEFAULT_MEAN_MQ_THR = 10
DEFAULT_MINIMUM_SEQUENCE_SIZE = 5980    # 20 % of the genome, 29903 bp


class Configuration:

    # file system storage
    ENV_COVIGATOR_STORAGE_FOLDER = "COVIGATOR_STORAGE_FOLDER"
    ENV_COVIGATOR_TEMP_FOLDER = "COVIGATOR_TEMP_FOLDER"
    ENV_COVIGATOR_DOWNLOAD_CONTENT_FOLDER = "COVIGATOR_DOWNLOAD_CONTENT_FOLDER"
    # database
    ENV_COVIGATOR_DB_PORT = "COVIGATOR_DB_PORT"
    ENV_COVIGATOR_DB_PASSWORD = "COVIGATOR_DB_PASSWORD"
    ENV_COVIGATOR_DB_USER = "COVIGATOR_DB_USER"
    ENV_COVIGATOR_DB_NAME = "COVIGATOR_DB_NAME"
    ENV_COVIGATOR_DB_HOST = "COVIGATOR_DB_HOST"
    ENV_COVIGATOR_DB_POOL_SIZE = "COVIGATOR_DB_POOL_SIZE"
    ENV_COVIGATOR_DB_MAX_OVERFLOW = "COVIGATOR_DB_MAX_OVERFLOW"
    ENV_COVIGATOR_TABLE_VERSION = "COVIGATOR_TABLE_VERSION"
    # dashboard
    ENV_COVIGATOR_DASHBOARD_HOST = "COVIGATOR_DASHBOARD_HOST"
    ENV_COVIGATOR_DASHBOARD_PORT = "COVIGATOR_DASHBOARD_PORT"
    # logs
    ENV_COVIGATOR_DASHBOARD_LOG_FILE = "COVIGATOR_DASHBOARD_LOG_FILE"
    ENV_COVIGATOR_PROCESSOR_LOG_FILE = "COVIGATOR_PROCESSOR_LOG_FILE"
    ENV_COVIGATOR_ACCESSOR_LOG_FILE = "COVIGATOR_ACCESSOR_LOG_FILE"
    # pipeline
    ENV_COVIGATOR_NEXTFLOW = "COVIGATOR_NEXTFLOW"
    ENV_COVIGATOR_WORKFLOW = "COVIGATOR_WORKFLOW"
    ENV_COVIGATOR_FORCE_PIPELINE = "COVIGATOR_FORCE_PIPELINE"
    ENV_COVIGATOR_SKIP_VCF_LOADING = "COVIGATOR_SKIP_VCF_LOADING"
    ENV_COVIGATOR_WORKFLOW_CPUS = "COVIGATOR_WORKFLOW_CPUS"
    ENV_COVIGATOR_WORKFLOW_MEMORY = "COVIGATOR_WORKFLOW_MEMORY"
    ENV_COVIGATOR_BATCH_SIZE = "COVIGATOR_BATCH_SIZE"
    ENV_COVIGATOR_RETRIES_DOWNLOAD = "COVIGATOR_RETRIES_DOWNLOAD"
    ENV_COVIGATOR_LOW_COVERAGE_THR = "COVIGATOR_LOW_COVERAGE_THR"
    ENV_COVIGATOR_SUBCLONAL_THR = "COVIGATOR_SUBCLONAL_THR"
    # references
    ENV_COVIGATOR_REF_FASTA = "COVIGATOR_REF_FASTA"
    # dask
    ENV_COVIGATOR_DASK_PORT = "COVIGATOR_DASK_PORT"
    # sample exclusion
    ENV_COVIGATOR_MEAN_MQ_THR = "COVIGATOR_MEAN_MQ_THR"
    ENV_COVIGATOR_MEAN_BQ_THR = "COVIGATOR_MEAN_BQ_THR"
    ENV_COVIGATOR_HORIZONTAL_COVERAGE_THR = "COVIGATOR_HORIZONTAL_COVERAGE_THR"
    ENV_COVIGATOR_MAX_SNVS = "COVIGATOR_MAX_SNVS"
    ENV_COVIGATOR_MAX_INSERTIONS = "COVIGATOR_MAX_INSERTIONS"
    ENV_COVIGATOR_MAX_DELETIONS = "COVIGATOR_MAX_DELETIONS"
    ENV_COVIGATOR_MIN_SEQUENCE_SIZE = "COVIGATOR_MIN_SEQUENCE_SIZE"

    def __init__(self, verbose=True):
        # local storage
        self.storage_folder = os.getenv(self.ENV_COVIGATOR_STORAGE_FOLDER, "/data/covigator")
        self.content_folder = os.getenv(self.ENV_COVIGATOR_DOWNLOAD_CONTENT_FOLDER)

        # database
        self.db_host = os.getenv(self.ENV_COVIGATOR_DB_HOST, "localhost")
        self.db_name = os.getenv(self.ENV_COVIGATOR_DB_NAME, "covigator")
        self.db_user = os.getenv(self.ENV_COVIGATOR_DB_USER, "covigator")
        self.db_password = os.getenv(self.ENV_COVIGATOR_DB_PASSWORD, "covigator")
        self.db_port = os.getenv(self.ENV_COVIGATOR_DB_PORT, "5432")
        self.db_pool_size = int(os.getenv(self.ENV_COVIGATOR_DB_POOL_SIZE, 5))
        self.db_max_overflow = int(os.getenv(self.ENV_COVIGATOR_DB_MAX_OVERFLOW, 10))
        self.db_table_version = os.environ.get(self.ENV_COVIGATOR_TABLE_VERSION, "_test")
        self.force_pipeline = os.environ.get(self.ENV_COVIGATOR_FORCE_PIPELINE, False)
        # this is useful for reloading metadata into the DB
        self.skip_vcf_loading = os.environ.get(self.ENV_COVIGATOR_SKIP_VCF_LOADING, False)

        # dashboard
        self.dash_host = os.getenv(self.ENV_COVIGATOR_DASHBOARD_HOST, "0.0.0.0")
        try:
            self.dash_port = int(os.getenv(self.ENV_COVIGATOR_DASHBOARD_PORT, "8050"))
        except ValueError as e:
            raise CovigatorDashBoardInitialisationError("The port needs to be a numeric value. " + str(e))

        # dask
        self.dask_port = os.getenv(self.ENV_COVIGATOR_DASK_PORT, "50218")

        # logs
        self.logfile_dash = os.getenv(self.ENV_COVIGATOR_DASHBOARD_LOG_FILE)
        self.logfile_accesor = os.getenv(self.ENV_COVIGATOR_ACCESSOR_LOG_FILE)
        self.logfile_processor = os.getenv(self.ENV_COVIGATOR_PROCESSOR_LOG_FILE)

        # references
        self.reference_genome = os.getenv(self.ENV_COVIGATOR_REF_FASTA)

        # pipeline
        self.nextflow = os.getenv(self.ENV_COVIGATOR_NEXTFLOW, "nextflow")
        self.workflow = os.getenv(self.ENV_COVIGATOR_WORKFLOW,
                                  "tron-bioinformatics/covigator-ngs-pipeline -r {version}".format(
                                      version=covigator.ANALYSIS_PIPELINE_VERSION))
        self.workflow_cpus = os.getenv(self.ENV_COVIGATOR_WORKFLOW_CPUS, "1")
        self.workflow_memory = os.getenv(self.ENV_COVIGATOR_WORKFLOW_MEMORY, "3g")
        self.batch_size = self.load_numeric_value(variable=self.ENV_COVIGATOR_BATCH_SIZE, default=1000)
        self.retries_download = self.load_numeric_value(variable=self.ENV_COVIGATOR_RETRIES_DOWNLOAD, default=3)
        self.low_coverage_threshold = self.load_float_value(
            variable=self.ENV_COVIGATOR_LOW_COVERAGE_THR, default=DEFAULT_LOW_COVERAGE_AF_THR)
        self.subclonal_threshold = self.load_float_value(
            variable=self.ENV_COVIGATOR_SUBCLONAL_THR, default=DEFAULT_SUBCLONAL_AF_THR)

        ## sample exclusion
        self.mean_mq_thr = self.load_numeric_value(variable=self.ENV_COVIGATOR_MEAN_MQ_THR, default=DEFAULT_MEAN_MQ_THR)
        self.mean_bq_thr = self.load_numeric_value(variable=self.ENV_COVIGATOR_MEAN_BQ_THR, default=DEFAULT_MEAN_BQ_THR)
        self.horizontal_coverage_thr = self.load_numeric_value(
            variable=self.ENV_COVIGATOR_HORIZONTAL_COVERAGE_THR, default=DEFAULT_HORIZONTAL_COVERAGE_THR)
        self.max_snvs = self.load_numeric_value(variable=self.ENV_COVIGATOR_MAX_SNVS, default=76)
        self.max_insertions = self.load_numeric_value(variable=self.ENV_COVIGATOR_MAX_INSERTIONS, default=10)
        self.max_deletions = self.load_numeric_value(variable=self.ENV_COVIGATOR_MAX_DELETIONS, default=10)
        self.min_sequence_size = self.load_numeric_value(variable=self.ENV_COVIGATOR_MIN_SEQUENCE_SIZE,
                                                         default=DEFAULT_MINIMUM_SEQUENCE_SIZE)

        # NOTE: the defaults are already set in the workflow config
        self.temp_folder = os.getenv(self.ENV_COVIGATOR_TEMP_FOLDER, "data/covigator-tmp")

        if verbose:
            self.log_configuration()

    def load_numeric_value(self, variable, default):
        try:
            value = int(os.getenv(variable, default))
        except ValueError as e:
            raise CovigatorDashBoardInitialisationError("{} needs to be a numeric value : {}".format(variable, str(e)))
        return value

    def load_float_value(self, variable, default):
        try:
            value = float(os.getenv(variable, default))
        except ValueError as e:
            raise CovigatorDashBoardInitialisationError("{} needs to be a float value : {}".format(variable, str(e)))
        return value

    def log_configuration(self):
        logger.info("Configuration")
        for k, v in self.__dict__.items():
            logger.info("{}={}".format(k, v))


def initialise_logs(logfile, sample_id: str = None):
    if logfile is not None:
        logzero.logfile(logfile)
    logzero.loglevel(logging.INFO)
    if sample_id is not None:
        logzero.formatter(logging.Formatter(
            '[%(levelname)1.1s %(asctime)s %(module)s:%(lineno)d ' + sample_id + '] %(message)s'))
    else:
        logzero.formatter(logging.Formatter(
            '[%(levelname)1.1s %(asctime)s %(module)s:%(lineno)d] %(message)s'))
