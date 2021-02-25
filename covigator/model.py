from datetime import datetime
from typing import List

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Float, Enum, DateTime, Integer, Boolean, Date, ForeignKey, UniqueConstraint, \
    ForeignKeyConstraint, BigInteger, JSON
import enum

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
    first_created = Column(Date)
    collection_date = Column(Date)
    instrument_platform = Column(String)
    instrument_model = Column(String)
    library_name = Column(String)
    nominal_length = Column(Integer)
    library_layout = Column(String)
    library_strategy = Column(String)
    library_source = Column(String)
    library_selection = Column(String)
    read_count = Column(BigInteger)
    base_count = Column(BigInteger)
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
    country_raw = Column(String)                                    # 'Denmark',
    country = Column(String)
    country_alpha_2 = Column(String)
    country_alpha_3 = Column(String)
    continent = Column(String)
    continent_alpha_2 = Column(String)

    def get_fastqs_ftp(self) -> List:
        return self.fastq_ftp.split(SEPARATOR) if self.fastq_ftp is not None else []

    def get_fastqs_md5(self) -> List:
        return self.fastq_md5.split(SEPARATOR) if self.fastq_md5 is not None else []


class Job(Base):

    __tablename__ = 'job'

    run_accession = Column(ForeignKey("ena_run.run_accession"), primary_key=True)

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
            fastq1 = next(filter(lambda x: "_1.fastq" in x, fastqs))
            fastq2 = next(filter(lambda x: "_2.fastq" in x, fastqs))
        else:
            fastq1 = fastqs[0]
            fastq2 = None
        return fastq1, fastq2


class Variant(Base):

    __tablename__ = 'variant'

    chromosome = Column(String, primary_key=True)
    position = Column(Integer, primary_key=True)
    reference = Column(String, primary_key=True)
    alternate = Column(String, primary_key=True)
    overlaps_multiple_genes = Column(Boolean, default=False)
    """
    ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
    Allele |                    C 
    Annotation |                missense_variant
    Annotation_Impact |         MODERATE 
    Gene_Name |                 orf1ab 
    Gene_ID |                   gene-GU280_gp01 
    Feature_Type |              transcript 
    Feature_ID |                TRANSCRIPT_gene-GU280_gp01 
    Transcript_BioType |        protein_coding 
    Rank |                      1/1 
    HGVS.c |                    c.33A>C 
    HGVS.p |                    p.Lys11Asn 
    cDNA.pos / cDNA.length |    33/21290 
    CDS.pos / CDS.length |      33/21290 
    AA.pos / AA.length |        11/7095 
    Distance | 
    ERRORS / WARNINGS / INFO'   
    """
    annotation = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)


class VariantObservation(Base):

    __tablename__ = 'variant_observation'

    sample = Column(ForeignKey("ena_run.run_accession"), primary_key=True)
    chromosome = Column(String, primary_key=True)
    position = Column(Integer, primary_key=True)
    reference = Column(String, primary_key=True)
    alternate = Column(String, primary_key=True)
    quality = Column(Float)
    filter = Column(String)
    """
    ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    ##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
    ##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
    ##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
    ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
    ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
    ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
    ##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
    ##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
    ##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
    ##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
    """
    dp = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    mapping_quality = Column(Integer)
    mapping_quality_zero_fraction = Column(Float)
    ForeignKeyConstraint(
        [chromosome, position, reference, alternate],
        [Variant.chromosome, Variant.position, Variant.reference, Variant.alternate])


class VariantCooccurrence(Base):

    __tablename__ = 'variant_cooccurrence'

    chromosome_one = Column(String, primary_key=True)
    position_one = Column(Integer, primary_key=True)
    reference_one = Column(String, primary_key=True)
    alternate_one = Column(String, primary_key=True)
    chromosome_two = Column(String, primary_key=True)
    position_two = Column(Integer, primary_key=True)
    reference_two = Column(String, primary_key=True)
    alternate_two = Column(String, primary_key=True)
    count = Column(Integer, default=0)

    ForeignKeyConstraint(
        [chromosome_one, position_one, reference_one, alternate_one],
        [Variant.chromosome, Variant.position, Variant.reference, Variant.alternate])
    ForeignKeyConstraint(
        [chromosome_two, position_two, reference_two, alternate_two],
        [Variant.chromosome, Variant.position, Variant.reference, Variant.alternate])


class Gene(Base):

    __tablename__ = 'gene'
    identifier = Column(String, primary_key=True)
    name = Column(String)
    data = Column(JSON)

