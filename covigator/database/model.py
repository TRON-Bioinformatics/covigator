import os
from datetime import datetime
from typing import List
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, String, Float, Enum, DateTime, Integer, Boolean, Date, ForeignKey, \
    ForeignKeyConstraint, BigInteger, JSON, Index
import enum

from covigator.configuration import Configuration


def get_table_versioned_name(basename, config: Configuration):
    return "{}{}".format(basename, config.db_table_version)


config = Configuration()
GENE_TABLE_NAME = get_table_versioned_name('gene', config=config)
LOG_TABLE_NAME = get_table_versioned_name('log', config=config)
VARIANT_COOCCURRENCE_TABLE_NAME = get_table_versioned_name('variant_cooccurrence', config=config)
VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('variant_observation', config=config)
VARIANT_TABLE_NAME = get_table_versioned_name('variant', config=config)
JOB_ENA_TABLE_NAME = get_table_versioned_name('job_ena', config=config)
JOB_GISAID_TABLE_NAME = get_table_versioned_name('job_gisaid', config=config)
SAMPLE_TABLE_NAME = get_table_versioned_name('sample', config=config)
SAMPLE_GISAID_TABLE_NAME = get_table_versioned_name('sample_gisaid', config=config)
SAMPLE_ENA_TABLE_NAME = get_table_versioned_name('sample_ena', config=config)
CONSERVATION_TABLE_NAME = get_table_versioned_name('conservation', config=config)
JOB_STATUS_CONSTRAINT_NAME = get_table_versioned_name('job_status', config=config)
DATA_SOURCE_CONSTRAINT_NAME = get_table_versioned_name('data_source', config=config)
COVIGATOR_MODULE_CONSTRAINT_NAME = get_table_versioned_name('covigator_module', config=config)
SEPARATOR = ";"

Base = declarative_base()


class Gene(Base):
    """
    This table holds the genes in the organism genome and its annotations in the field `data`. The annotations
    are in JSON format. This data is fetched from ftp://ftp.ensemblgenomes.org/pub/viruses/json/sars_cov_2/sars_cov_2.json
    """
    __tablename__ = GENE_TABLE_NAME

    identifier = Column(String, primary_key=True)
    name = Column(String)
    start = Column(Integer, index=True)
    end = Column(Integer)
    data = Column(JSON)

    def get_pfam_domains(self):
        protein_features = self.data.get("transcripts", [])[0].get("translations", [])[0].get("protein_features")
        pfam_protein_features = [f for f in protein_features if f.get("dbname") == "Pfam"]
        return sorted(pfam_protein_features, key=lambda d: int(d.get("start")))


class JobStatus(enum.Enum):
    """
    Valid job status
    """
    __constraint_name__ = JOB_STATUS_CONSTRAINT_NAME

    PENDING = 1
    QUEUED = 2
    DOWNLOADED = 3
    PROCESSED = 4
    LOADED = 5
    FAILED_DOWNLOAD = 6
    FAILED_PROCESSING = 7
    FAILED_LOAD = 8
    COOCCURRENCE = 9
    FAILED_COOCCURRENCE = 10
    FINISHED = 11


class DataSource(enum.Enum):
    """
    Valid sources of data
    """
    __constraint_name__ = DATA_SOURCE_CONSTRAINT_NAME

    ENA = 1
    GISAID = 2


class SampleGisaid(Base):
    """
    The table that holds all metadata for a GISAID sample
    """
    __tablename__ = SAMPLE_GISAID_TABLE_NAME

    run_accession = Column(String, primary_key=True)
    date = Column(Date)
    # Host information
    host_tax_id = Column(String)
    host = Column(String)
    # geographical data
    country_raw = Column(String)
    region = Column(String)
    country = Column(String)
    country_alpha_2 = Column(String)
    country_alpha_3 = Column(String)
    continent = Column(String)
    continent_alpha_2 = Column(String)
    site = Column(String)
    site2 = Column(String)
    sequence = Column(JSON)


class SampleEna(Base):
    """
    The table that holds all metadata for a ENA sample
    """
    __tablename__ = SAMPLE_ENA_TABLE_NAME

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


class Sample(Base):
    """
    This table holds all samples loaded into Covigator irrespective of the data source.
    The same sample may be loaded from different data sources.
    There are foreign keys fields pointing to the source-specific tables with all metadata for the sample.
    """
    __tablename__ = SAMPLE_TABLE_NAME

    id = Column(String, primary_key=True)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__), primary_key=True)
    # NOTE: should have only one filled, either ena_id or gisaid_id and be coherent with the value of source
    ena_id = Column(ForeignKey("{}.run_accession".format(SampleEna.__tablename__)))
    gisaid_id = Column(ForeignKey("{}.run_accession".format(SampleGisaid.__tablename__)))


class JobGisaid(Base):
    """
    The table that holds an GISAID job
    """
    __tablename__ = JOB_GISAID_TABLE_NAME

    run_accession = Column(ForeignKey("{}.run_accession".format(SampleGisaid.__tablename__)), primary_key=True)

    # job status
    status = Column(Enum(JobStatus), default=JobStatus.PENDING)
    created_at = Column(DateTime(timezone=True), nullable=False, default=datetime.now())
    queued_at = Column(DateTime(timezone=True))
    #downloaded_at = Column(DateTime(timezone=True))
    analysed_at = Column(DateTime(timezone=True))
    #cleaned_at = Column(DateTime(timezone=True))
    loaded_at = Column(DateTime(timezone=True))
    failed_at = Column(DateTime(timezone=True))
    error_message = Column(String)

    # TODO: Test DB size, sequence should be approx 30000bp, compress maybe
    # local files storage
    vcf_path = Column(String)


class JobEna(Base):
    """
    The table that holds an ENA job
    """
    __tablename__ = JOB_ENA_TABLE_NAME

    run_accession = Column(ForeignKey("{}.run_accession".format(SampleEna.__tablename__)), primary_key=True)

    # job status
    status = Column(Enum(JobStatus, name=JobStatus.__constraint_name__), default=JobStatus.PENDING)
    created_at = Column(DateTime(timezone=True), nullable=False, default=datetime.now())
    queued_at = Column(DateTime(timezone=True))
    downloaded_at = Column(DateTime(timezone=True))
    analysed_at = Column(DateTime(timezone=True))
    cleaned_at = Column(DateTime(timezone=True))
    loaded_at = Column(DateTime(timezone=True))
    cooccurrence_at = Column(DateTime(timezone=True))
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
    """
    A variant with its specific annotations. THis does not contain any sample specific annotations.
    """
    __tablename__ = VARIANT_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
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

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class VariantObservation(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = VARIANT_OBSERVATION_TABLE_NAME

    source = Column(Enum(DataSource, name=DataSource.__constraint_name__), primary_key=True)
    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
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
    ForeignKeyConstraint([sample, source], [Sample.id, Sample.source])
    ForeignKeyConstraint([variant_id], [Variant.variant_id])


class VariantCooccurrence(Base):

    __tablename__ = VARIANT_COOCCURRENCE_TABLE_NAME

    variant_id_one = Column(String, primary_key=True)
    variant_id_two = Column(String, primary_key=True)
    count = Column(Integer, default=0)

    ForeignKeyConstraint([variant_id_one], [Variant.variant_id])
    ForeignKeyConstraint([variant_id_two], [Variant.variant_id])


class CovigatorModule(enum.Enum):

    __constraint_name__ = COVIGATOR_MODULE_CONSTRAINT_NAME

    ACCESSOR = 1
    PROCESSOR = 2


class Log(Base):
    """
    The table that holds an ENA job
    """
    __tablename__ = LOG_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)

    start = Column(DateTime(timezone=True), nullable=False)
    end = Column(DateTime(timezone=True), nullable=False)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__), nullable=False)
    module = Column(Enum(CovigatorModule, name=CovigatorModule.__constraint_name__), nullable=False)

    has_error = Column(Boolean, default=False)
    error_message = Column(String)
    processed = Column(Integer)
    data = Column(JSON)


class Conservation(Base):
    """
    Table to hold conservation scores
    """
    __tablename__ = CONSERVATION_TABLE_NAME

    chromosome = Column(String, primary_key=True)
    start = Column(Integer, primary_key=True)
    end = Column(Integer, primary_key=True)
    conservation = Column(Float)
    conservation_sarbecovirus = Column(Float)
    conservation_vertebrates = Column(Float)
