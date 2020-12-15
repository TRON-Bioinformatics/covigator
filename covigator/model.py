from typing import List

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Float, Enum
import enum


Base = declarative_base()


class EnaRunStatus(enum.Enum):
    NEW = 1
    DOWNLOADING = 2
    DOWNLOADED = 3
    FAILED = 4


class EnaRun(Base):

    __tablename__ = 'ena_run'

    # data on run
    run_accession = Column(String, primary_key=True)            # 'ERR4080473',
    sample_accession = Column(String)                           # 'SAMEA6798401',
    scientific_name = Column(String)
    study_accession = Column(String)
    experiment_accession = Column(String)
    first_created = Column(String)
    collection_date = Column(String)
    instrument_platform = Column(String)
    instrument_model = Column(String)
    sample_collection = Column(String)
    sequencing_method = Column(String)
    center_name = Column(String)

    # FASTQs
    fastq_ftp = Column(String, nullable=False)                  # 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080473/ERR4080473_1.fastq.gz',
    fastq_md5 = Column(String, nullable=False)

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

    # covigator status
    status = Column(Enum(EnaRunStatus), default=EnaRunStatus.NEW)

    def get_fastqs_ftp(self) -> List:
        return self.fastq_ftp.split(";")
