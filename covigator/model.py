from datetime import datetime
from typing import List
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Float, Enum, create_engine, DateTime, Integer
from sqlalchemy.orm import Session, sessionmaker
import enum
import os
from logzero import logger

from covigator import ENV_COVIGATOR_DB_PORT, ENV_COVIGATOR_DB_PASSWORD, ENV_COVIGATOR_DB_USER, ENV_COVIGATOR_DB_NAME, \
    ENV_COVIGATOR_DB_HOST

SEPARATOR = ";"

Base = declarative_base()


class JobStatus(enum.Enum):
    PENDING = 1
    QUEUED = 2
    DOWNLOADED = 3
    PROCESSED = 4
    LOADED = 5
    FAILED_DOWNLOAD = 6
    FAILED_PROCESSING = 7
    FAILED_LOAD = 8


class EnaRun(Base):

    __tablename__ = 'ena_run'

    # data on run
    # TODO: add foreign keys to jobs
    run_accession = Column(String, primary_key=True)            # 'ERR4080473',
    sample_accession = Column(String)                           # 'SAMEA6798401',
    scientific_name = Column(String)
    study_accession = Column(String)
    experiment_accession = Column(String)
    # TODO: store these as dates
    first_created = Column(String)
    collection_date = Column(String)
    instrument_platform = Column(String)
    instrument_model = Column(String)
    library_name = Column(String)
    nominal_length = Column(String)
    library_layout = Column(String)
    library_strategy = Column(String)
    library_source = Column(String)
    library_selection = Column(String)
    read_count = Column(String)
    base_count = Column(String)
    sample_collection = Column(String)
    sequencing_method = Column(String)
    center_name = Column(String)

    # FASTQs
    fastq_ftp = Column(String, nullable=False)                  # 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080473/ERR4080473_1.fastq.gz',
    fastq_md5 = Column(String, nullable=False)
    num_fastqs = Column(Integer, nullable=False)

    # data on host
    host_tax_id = Column(String)                                # '9606',
    host_sex = Column(String)
    host_body_site = Column(String)
    host_gravidity = Column(String)
    host_phenotype = Column(String)
    host_genotype = Column(String)

    # geographical data
    lat = Column(Float)
    lon = Column(Float)
    country = Column(String)                                    # 'Denmark',

    def get_fastqs_ftp(self) -> List:
        return self.fastq_ftp.split(SEPARATOR) if self.fastq_ftp is not None else []

    def get_fastqs_md5(self) -> List:
        return self.fastq_md5.split(SEPARATOR) if self.fastq_md5 is not None else []


class Job(Base):

    __tablename__ = 'job'

    run_accession = Column(String, primary_key=True)

    # job status
    status = Column(Enum(JobStatus), default=JobStatus.PENDING)
    created_at = Column(DateTime(timezone=True), nullable=False, default=datetime.now())
    queued_at = Column(DateTime(timezone=True))
    downloaded_at = Column(DateTime(timezone=True))
    analysed_at = Column(DateTime(timezone=True))
    cleaned_at = Column(DateTime(timezone=True))
    loaded_at = Column(DateTime(timezone=True))
    failed_at = Column(DateTime(timezone=True))
    error_message = Column(String)

    # local files storage
    fastq_path = Column(String)     # the local path where FASTQ files are stored in semi colon separated list
    vcf_path = Column(String)

    def get_fastq_paths(self):
        return self.fastq_path.split(SEPARATOR) if self.fastq_path is not None else []

    def get_fastq1_and_fastq2(self):
        fastqs = self.get_fastq_paths()
        if len(fastqs) > 1:
            fastq1 = next("_1.fastq" in x for x in fastqs)
            fastq2 = next("_2.fastq" in x for x in fastqs)
        else:
            fastq1 = fastqs[0]
            fastq2 = None
        return fastq1, fastq2


class Database:

    def __init__(self, test=False):
        if test:
            # creates a SQLite in memory database for testing purposes
            self.engine = create_engine('sqlite://')
        else:
            host = os.getenv(ENV_COVIGATOR_DB_HOST, "0.0.0.0")
            database = os.getenv(ENV_COVIGATOR_DB_NAME, "covigator")
            user = os.getenv(ENV_COVIGATOR_DB_USER, "covigator")
            password = os.getenv(ENV_COVIGATOR_DB_PASSWORD, "covigator")
            port = os.getenv(ENV_COVIGATOR_DB_PORT, "5432")
            db_uri = "postgresql+psycopg2://%s:%s@%s:%s/%s" % (
                user,
                password,
                host,
                port,
                database,
            )
            self.engine = create_engine(db_uri)
        self.engine.connect()
        self.Session = sessionmaker(bind=self.engine, autoflush=False)
        self.create_database()

    def create_database(self):
        # this creates all tables in the database (when it exists nothing happens)
        Base.metadata.create_all(self.engine)
        logger.info("Database initialized")

    def get_database_session(self) -> Session:
        return self.Session()
