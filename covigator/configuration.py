import logging
import os
import logzero
from covigator.exceptions import CovigatorDashBoardInitialisationError


class Configuration:

    # configuration environment variables
    ENV_COVIGATOR_STORAGE_FOLDER = "COVIGATOR_STORAGE_FOLDER"
    ENV_COVIGATOR_DB_PORT = "COVIGATOR_DB_PORT"
    ENV_COVIGATOR_DB_PASSWORD = "COVIGATOR_DB_PASSWORD"
    ENV_COVIGATOR_DB_USER = "COVIGATOR_DB_USER"
    ENV_COVIGATOR_DB_NAME = "COVIGATOR_DB_NAME"
    ENV_COVIGATOR_DB_HOST = "COVIGATOR_DB_HOST"
    ENV_COVIGATOR_DB_POOL_SIZE = "COVIGATOR_DB_POOL_SIZE"
    ENV_COVIGATOR_DB_MAX_OVERFLOW = "COVIGATOR_DB_MAX_OVERFLOW"
    ENV_COVIGATOR_DASHBOARD_HOST = "COVIGATOR_DASHBOARD_HOST"
    ENV_COVIGATOR_DASHBOARD_PORT = "COVIGATOR_DASHBOARD_PORT"
    ENV_COVIGATOR_DASHBOARD_LOG_FILE = "COVIGATOR_DASHBOARD_LOG_FILE"
    ENV_COVIGATOR_PROCESSOR_LOG_FILE = "COVIGATOR_PROCESSOR_LOG_FILE"
    ENV_COVIGATOR_ACCESSOR_LOG_FILE = "COVIGATOR_ACCESSOR_LOG_FILE"
    ENV_COVIGATOR_BIN_FASTP = "COVIGATOR_BIN_FASTP"
    ENV_COVIGATOR_BIN_SAMTOOLS = "COVIGATOR_BIN_SAMTOOLS"
    ENV_COVIGATOR_BIN_BWA = "COVIGATOR_BIN_BWA"
    ENV_COVIGATOR_BIN_BCFTOOLS = "COVIGATOR_BIN_BCFTOOLS"
    ENV_COVIGATOR_BIN_SNPEFF = "COVIGATOR_BIN_SNPEFF"
    ENV_COVIGATOR_BIN_BGZIP = "COVIGATOR_BIN_BGZIP"
    ENV_COVIGATOR_BIN_TABIX = "COVIGATOR_BIN_TABIX"
    ENV_COVIGATOR_BIN_JAVA = "COVIGATOR_BIN_JAVA"
    ENV_COVIGATOR_REF_FASTA = "COVIGATOR_REF_FASTA"
    ENV_COVIGATOR_REF_GISAID = "COVIGATOR_REF_GISAID"
    ENV_COVIGATOR_SEQ_GISAID = "COVIGATOR_SEQ_GISAID"
    ENV_COVIGATOR_META_GISAID = "COVIGATOR_META_GISAID"
    ENV_COVIGATOR_DASK_PORT = "COVIGATOR_DASK_PORT"

    ENV_COVIGATOR_TABLE_VERSION = "COVIGATOR_TABLE_VERSION"

    def __init__(self):
        # local storage
        self.storage_folder = os.getenv(self.ENV_COVIGATOR_STORAGE_FOLDER, "./data/covigator")

        # database
        self.db_host = os.getenv(self.ENV_COVIGATOR_DB_HOST, "0.0.0.0")
        self.db_name = os.getenv(self.ENV_COVIGATOR_DB_NAME, "covigator")
        self.db_user = os.getenv(self.ENV_COVIGATOR_DB_USER, "covigator")
        self.db_password = os.getenv(self.ENV_COVIGATOR_DB_PASSWORD, "covigator")
        self.db_port = os.getenv(self.ENV_COVIGATOR_DB_PORT, "5432")
        self.db_pool_size = int(os.getenv(self.ENV_COVIGATOR_DB_POOL_SIZE, 5))
        self.db_max_overflow = int(os.getenv(self.ENV_COVIGATOR_DB_MAX_OVERFLOW, 10))
        self.db_table_version = os.environ.get(self.ENV_COVIGATOR_TABLE_VERSION, "")

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

        # pipeline binaries
        self.fastp = os.getenv(self.ENV_COVIGATOR_BIN_FASTP, "/code/fastp/fastp")
        self.samtools = os.getenv(self.ENV_COVIGATOR_BIN_SAMTOOLS, "/code/samtools/1.9/samtools")
        self.bwa = os.getenv(self.ENV_COVIGATOR_BIN_BWA, "/code/bwa/0.7.17/bwa")
        self.bcftools = os.getenv(self.ENV_COVIGATOR_BIN_BCFTOOLS, "/code/bcftools/1.9/bcftools")
        self.snpeff = os.getenv(self.ENV_COVIGATOR_BIN_SNPEFF, "/code/snpEFF_latest_core/snpEff/snpEff.jar")
        self.bgzip = os.getenv(self.ENV_COVIGATOR_BIN_BGZIP, "bgzip")
        self.tabix = os.getenv(self.ENV_COVIGATOR_BIN_TABIX, "tabix")
        self.java = os.getenv(self.ENV_COVIGATOR_BIN_JAVA, "java")

        # references
        self.reference_genome = os.getenv(self.ENV_COVIGATOR_REF_FASTA,
                                          "/scratch/info/projects/SARS-CoV-2/index/MN908947.3.fa")


def initialise_logs(logfile, sample_id: str = None):
    if logfile is not None:
        logzero.logfile(logfile, maxBytes=1e6, backupCount=3)
    logzero.loglevel(logging.INFO)
    if sample_id is not None:
        logzero.formatter(logging.Formatter('%(color)s[%(levelname)1.1s %(asctime)s %(module)s:%(lineno)d ' + sample_id +
                          ']%(end_color)s %(message)s'))