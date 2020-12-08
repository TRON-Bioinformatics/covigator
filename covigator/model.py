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

    run_accession = Column(String, primary_key=True)            # 'ERR4080473',
    sample_accession = Column(String)                           # 'SAMEA6798401',
    fastq_ftp = Column(String, nullable=False)                  # 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080473/ERR4080473_1.fastq.gz',
    host_tax_id = Column(String)                                # '9606',
    host_sex = Column(String)
    lat = Column(Float)
    lon = Column(Float)
    country = Column(String)                                    # 'Denmark',
    status = Column(Enum(EnaRunStatus), default=EnaRunStatus.NEW)
    file_path = Column(String)

    def get_fastqs_ftp(self) -> List:
        return self.fastq_ftp.split(";")
