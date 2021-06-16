import logging
import os
import logzero
from covigator.exceptions import CovigatorDashBoardInitialisationError


class Configuration:

    # file system storage
    ENV_COVIGATOR_STORAGE_FOLDER = "COVIGATOR_STORAGE_FOLDER"
    ENV_COVIGATOR_TEMP_FOLDER = "COVIGATOR_TEMP_FOLDER"
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
    ENV_COVIGATOR_TRONFLOW_BWA = "COVIGATOR_TRONFLOW_BWA"
    ENV_COVIGATOR_TRONFLOW_BAM_PREPROCESSING = "COVIGATOR_TRONFLOW_BAM_PREPROCESSING"
    ENV_COVIGATOR_TRONFLOW_VARIANT_NORMALIZATION = "COVIGATOR_TRONFLOW_VARIANT_NORMALIZATION"
    ENV_COVIGATOR_FORCE_PIPELINE = "COVIGATOR_FORCE_PIPELINE"

    # references
    ENV_COVIGATOR_REF_FASTA = "COVIGATOR_REF_FASTA"
    ENV_COVIGATOR_GENE_ANNOTATIONS = "COVIGATOR_GENE_ANNOTATIONS"
    # dask
    ENV_COVIGATOR_DASK_PORT = "COVIGATOR_DASK_PORT"

    def __init__(self):
        # local storage
        self.storage_folder = os.getenv(self.ENV_COVIGATOR_STORAGE_FOLDER, "/data/covigator")

        # database
        self.db_host = os.getenv(self.ENV_COVIGATOR_DB_HOST, "0.0.0.0")
        self.db_name = os.getenv(self.ENV_COVIGATOR_DB_NAME, "covigator")
        self.db_user = os.getenv(self.ENV_COVIGATOR_DB_USER, "covigator")
        self.db_password = os.getenv(self.ENV_COVIGATOR_DB_PASSWORD, "covigator")
        self.db_port = os.getenv(self.ENV_COVIGATOR_DB_PORT, "5432")
        self.db_pool_size = int(os.getenv(self.ENV_COVIGATOR_DB_POOL_SIZE, 5))
        self.db_max_overflow = int(os.getenv(self.ENV_COVIGATOR_DB_MAX_OVERFLOW, 10))
        self.db_table_version = os.environ.get(self.ENV_COVIGATOR_TABLE_VERSION, "")
        self.force_pipeline = os.environ.get(self.ENV_COVIGATOR_FORCE_PIPELINE, False)

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
        self.reference_gene_annotations = os.getenv(self.ENV_COVIGATOR_GENE_ANNOTATIONS)

        self.nextflow = os.getenv(self.ENV_COVIGATOR_NEXTFLOW, "nextflow")
        self.workflow = os.getenv(self.ENV_COVIGATOR_WORKFLOW, "tron-bioinformatics/covigator-ngs-pipeline -r v0.3.0")
        # NOTE: the defaults are already set in the workflow config
        self.tronflow_bwa = os.getenv(self.ENV_COVIGATOR_TRONFLOW_BWA)
        self.tronflow_bam_preprocessing = os.getenv(self.ENV_COVIGATOR_TRONFLOW_BAM_PREPROCESSING)
        self.tronflow_variant_normalization = os.getenv(self.ENV_COVIGATOR_TRONFLOW_VARIANT_NORMALIZATION)
        self.temp_folder = os.getenv(self.ENV_COVIGATOR_TEMP_FOLDER, "/data/covigator-tmp")


def initialise_logs(logfile, sample_id: str = None):
    if logfile is not None:
        logzero.logfile(logfile)
    logzero.loglevel(logging.INFO)
    if sample_id is not None:
        logzero.formatter(logging.Formatter(
            '[%(levelname)1.1s %(asctime)s %(module)s:%(lineno)d ' + sample_id + '] %(message)s'))
