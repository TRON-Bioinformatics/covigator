import os
from datetime import datetime
from typing import List
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import Column, String, Float, Enum, DateTime, Integer, Boolean, Date, ForeignKey, \
    ForeignKeyConstraint, BigInteger, JSON, Index
import enum
from covigator.configuration import Configuration


def get_table_versioned_name(basename, config: Configuration):
    return "{}{}".format(basename, config.db_table_version)


config = Configuration()
GENE_TABLE_NAME = get_table_versioned_name('gene', config=config)
DOMAIN_TABLE_NAME = get_table_versioned_name('domain', config=config)
LOG_TABLE_NAME = get_table_versioned_name('log', config=config)
LAST_UPDATE_TABLE_NAME = get_table_versioned_name('last_update', config=config)
VARIANT_COOCCURRENCE_TABLE_NAME = get_table_versioned_name('variant_cooccurrence', config=config)
COVID19_PORTAL_VARIANT_COOCCURRENCE_TABLE_NAME = get_table_versioned_name('covid19portal_variant_cooccurrence', config=config)
VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('variant_observation', config=config)
SUBCLONAL_VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('subclonal_variant_observation', config=config)
LOW_FREQUENCY_VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('low_frequency_variant_observation', config=config)
LOW_QUALITY_CLONAL_VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('lq_clonal_variant_observation', config=config)
COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME = get_table_versioned_name('variant_observation_covid19portal', config=config)
VARIANT_TABLE_NAME = get_table_versioned_name('variant', config=config)
SUBCLONAL_VARIANT_TABLE_NAME = get_table_versioned_name('subclonal_variant', config=config)
LOW_QUALITY_CLONAL_VARIANT_TABLE_NAME = get_table_versioned_name('lq_clonal_variant', config=config)
LOW_FREQUENCY_VARIANT_TABLE_NAME = get_table_versioned_name('low_frequency_variant', config=config)
COVID19_PORTAL_VARIANT_TABLE_NAME = get_table_versioned_name('variant_covid19portal', config=config)
SAMPLE_ENA_TABLE_NAME = get_table_versioned_name('sample_ena', config=config)
SAMPLE_COVID19_PORTAL_TABLE_NAME = get_table_versioned_name('sample_covid19_portal', config=config)
CONSERVATION_TABLE_NAME = get_table_versioned_name('conservation', config=config)
PRECOMPUTED_VARIANTS_PER_SAMPLE_TABLE_NAME = get_table_versioned_name('precomputed_variants_per_sample', config=config)
PRECOMPUTED_SUBSTITUTIONS_COUNTS_TABLE_NAME = get_table_versioned_name('precomputed_substitutions_counts', config=config)
PRECOMPUTED_INDEL_LENGTH_TABLE_NAME = get_table_versioned_name('precomputed_indel_length', config=config)
PRECOMPUTED_ANNOTATION_TABLE_NAME = get_table_versioned_name('precomputed_annotation', config=config)
PRECOMPUTED_OCCURRENCE_TABLE_NAME = get_table_versioned_name('precomputed_top_occurrence', config=config)
PRECOMPUTED_NS_S_COUNTS_TABLE_NAME = get_table_versioned_name('precomputed_ns_s_counts', config=config)
PRECOMPUTED_DN_DS_BY_DOMAIN_TABLE_NAME = get_table_versioned_name('precomputed_dn_ds_by_domain', config=config)
PRECOMPUTED_TABLE_COUNTS_TABLE_NAME = get_table_versioned_name('precomputed_table_counts', config=config)
PRECOMPUTED_VARIANT_ABUNDANCE_HIST_TABLE_NAME = get_table_versioned_name('precomputed_variant_abundance_histogram', config=config)
PRECOMPUTED_VARIANTS_PER_LINEAGE_TABLE_NAME = get_table_versioned_name('precomputed_variants_per_lineage', config=config)
JOB_STATUS_CONSTRAINT_NAME = get_table_versioned_name('job_status', config=config)
DATA_SOURCE_CONSTRAINT_NAME = get_table_versioned_name('data_source', config=config)
COVIGATOR_MODULE_CONSTRAINT_NAME = get_table_versioned_name('covigator_module', config=config)
REGION_TYPE_CONSTRAINT_NAME = get_table_versioned_name('region_type', config=config)
VARIANT_TYPE_CONSTRAINT_NAME = get_table_versioned_name('variant_type', config=config)
LINEAGE_TABLE_NAME = get_table_versioned_name('lineage', config=config)
CONSTELLATION_SITES_TABLE_NAME = get_table_versioned_name('lineage_defining_variant', config=config)
LINEAGE_SITES_JUNCTION_TABLE_NAME = get_table_versioned_name('lineage_variant', config=config)
VARIANT_LEVEL_CONSTRAINT_NAME = get_table_versioned_name('variant_level', config=config)
NEWS_SECTION_TABLE_NAME = get_table_versioned_name('news_section', config=config)
NEWS_TYPE_CONSTRAINT_NAME = get_table_versioned_name('news_type', config=config)
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
    fraction_synonymous = Column(Float)
    fraction_non_synonymous = Column(Float)


class Domain(Base):
    """
    This table holds the Pfam domains for the genes
    """
    __tablename__ = DOMAIN_TABLE_NAME

    name = Column(String, primary_key=True)
    description = Column(String)
    start = Column(Integer, index=True)
    end = Column(Integer)
    fraction_synonymous = Column(Float)
    fraction_non_synonymous = Column(Float)
    gene_identifier = Column(ForeignKey("{}.identifier".format(Gene.__tablename__)))
    gene_name = Column(String)


class JobStatus(enum.Enum):
    """
    Valid job status
    """
    __constraint_name__ = JOB_STATUS_CONSTRAINT_NAME

    PENDING = 1
    FAILED_DOWNLOAD = 2
    DOWNLOADED = 3
    QUEUED = 4
    FAILED_PROCESSING = 5
    FINISHED = 6
    HOLD = 7
    EXCLUDED = 8


class DataSource(enum.Enum):
    """
    Valid sources of data
    """
    __constraint_name__ = DATA_SOURCE_CONSTRAINT_NAME

    ENA = 1
    GISAID = 2
    COVID19_PORTAL = 3


class VariantType(enum.Enum):
    __constraint_name__ = VARIANT_TYPE_CONSTRAINT_NAME

    SNV = 1
    INSERTION = 2
    DELETION = 3
    MNV = 4
    COMPLEX = 5


class SampleEna(Base):
    """
    The table that holds all metadata for a ENA sample
    """
    __tablename__ = SAMPLE_ENA_TABLE_NAME

    # data on run
    run_accession = Column(String, primary_key=True)            # 'ERR4080473',

    # DEPRECATED
    finished = Column(Boolean)
    ####

    sample_accession = Column(String)                           # 'SAMEA6798401',
    scientific_name = Column(String)
    study_accession = Column(String)
    experiment_accession = Column(String)
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

    # counts of variants
    count_snvs = Column(Integer)
    count_insertions = Column(Integer)
    count_deletions = Column(Integer)
    count_subclonal_snvs = Column(Integer)
    count_subclonal_insertions = Column(Integer)
    count_subclonal_deletions = Column(Integer)
    count_low_frequency_snvs = Column(Integer)
    count_low_frequency_insertions = Column(Integer)
    count_low_frequency_deletions = Column(Integer)

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
    sample_folder = Column(String)
    fastq_path = Column(String)  # the local path where FASTQ files are stored in semi colon separated list
    lofreq_vcf_path = Column(String)
    ivar_vcf_path = Column(String)
    gatk_vcf_path = Column(String)
    bcftools_vcf_path = Column(String)
    lofreq_pangolin_path = Column(String)
    ivar_pangolin_path = Column(String)
    gatk_pangolin_path = Column(String)
    bcftools_pangolin_path = Column(String)
    fastp_path = Column(String)
    deduplication_metrics_path = Column(String)     # DEPRECATED
    horizontal_coverage_path = Column(String)
    vertical_coverage_path = Column(String)

    # FASTP results
    qc = Column(JSON)

    # horizontal coverage analysis results
    num_reads = Column(Integer)
    covered_bases = Column(Integer)
    coverage = Column(Float)
    mean_depth = Column(Float)
    mean_base_quality = Column(Float)
    mean_mapping_quality = Column(Float)

    # pango output (corresponding only to LoFreq)
    pangolin_lineage = Column(String)
    pangolin_conflict = Column(Float)
    pangolin_ambiguity_score = Column(Float)
    pangolin_scorpio_call = Column(String)
    pangolin_scorpio_support = Column(Float)
    pangolin_scorpio_conflict = Column(Float)
    pangolin_version = Column(String)
    pangolin_pangolin_version = Column(String)
    pangolin_scorpio_version = Column(String)
    pangolin_constellation_version = Column(String)
    pangolin_qc_status = Column(String)
    pangolin_qc_notes = Column(String)
    pangolin_note = Column(String)

    # Picard deduplicatio output   ## DEPRECATED
    percent_duplication = Column(Float)
    unpaired_reads_examined = Column(Integer)
    read_pairs_examined = Column(Integer)
    secondary_or_supplementary_reads = Column(Integer)
    unmapped_reads = Column(Integer)
    unpaired_read_duplicates = Column(Integer)
    read_pair_duplicates = Column(Integer)
    read_pair_optical_duplicates = Column(Integer)

    # intrahost sample filters
    intrahost_filter = Column(Boolean, default=False)
    potential_coinfection = Column(Boolean, default=False)

    covigator_accessor_version = Column(String)
    covigator_processor_version = Column(String)

    def get_sample_folder(self, base_folder):
        return os.path.join(
            base_folder,
            self.collection_date.strftime("%Y%m%d") if self.collection_date is not None else "nodate",
            self.run_accession)

    def get_fastqs_ftp(self) -> List:
        return self.fastq_ftp.split(SEPARATOR) if self.fastq_ftp is not None else []

    def get_fastqs_md5(self) -> List:
        return self.fastq_md5.split(SEPARATOR) if self.fastq_md5 is not None else []

    def get_fastq_paths(self) -> List:
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


class SampleCovid19Portal(Base):
    """
    The table that holds all metadata for a ENA sample
    """
    __tablename__ = SAMPLE_COVID19_PORTAL_TABLE_NAME

    run_accession = Column(String, primary_key=True)
    first_created = Column(Date)
    collection_date = Column(Date)
    last_modification_date = Column(Date)
    center_name = Column(String)
    isolate = Column(String)
    molecule_type = Column(String)

    # geographical data
    country_raw = Column(String)  # 'Denmark',
    country = Column(String)
    country_alpha_2 = Column(String)
    country_alpha_3 = Column(String)
    continent = Column(String)
    continent_alpha_2 = Column(String)

    # FASTA
    fasta_url = Column(String, nullable=False)

    # local files storage
    sample_folder = Column(String)
    fasta_path = Column(String)
    vcf_path = Column(String)
    pangolin_path = Column(String)

    # sequence  information
    sequence_length = Column(Integer)
    count_n_bases = Column(Integer)
    count_ambiguous_bases = Column(Integer)

    # counts of variants
    count_snvs = Column(Integer)
    count_insertions = Column(Integer)
    count_deletions = Column(Integer)

    # pango output
    pangolin_lineage = Column(String)
    pangolin_conflict = Column(Float)
    pangolin_ambiguity_score = Column(Float)
    pangolin_scorpio_call = Column(String)
    pangolin_scorpio_support = Column(Float)
    pangolin_scorpio_conflict = Column(Float)
    pangolin_version = Column(String)
    pangolin_pangolin_version = Column(String)
    pangolin_scorpio_version = Column(String)
    pangolin_constellation_version = Column(String)
    pangolin_qc_status = Column(String)
    pangolin_qc_notes = Column(String)
    pangolin_note = Column(String)

    # job status
    status = Column(Enum(JobStatus, name=JobStatus.__constraint_name__), default=JobStatus.PENDING)
    created_at = Column(DateTime(timezone=True), nullable=False, default=datetime.now())
    queued_at = Column(DateTime(timezone=True))
    analysed_at = Column(DateTime(timezone=True))
    loaded_at = Column(DateTime(timezone=True))
    failed_at = Column(DateTime(timezone=True))
    error_message = Column(String)

    covigator_accessor_version = Column(String)
    covigator_processor_version = Column(String)

    def get_sample_folder(self, base_folder: str):
        return os.path.join(
            base_folder,
            self.collection_date.strftime("%Y%m%d") if self.collection_date is not None else "nodate",
            self.run_accession)

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
    SnpEff annotations
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
    annotation_highest_impact = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class VariantCovid19Portal(Base):
    """
    A variant with its specific annotations. THis does not contain any sample specific annotations.
    """
    __tablename__ = COVID19_PORTAL_VARIANT_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
    overlaps_multiple_genes = Column(Boolean, default=False)
    """
    SnpEff annotations
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
    annotation_highest_impact = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class SubclonalVariant(Base):

    __tablename__ = SUBCLONAL_VARIANT_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
    overlaps_multiple_genes = Column(Boolean, default=False)

    annotation = Column(String, index=True)
    annotation_highest_impact = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class LowFrequencyVariant(Base):

    __tablename__ = LOW_FREQUENCY_VARIANT_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
    overlaps_multiple_genes = Column(Boolean, default=False)

    annotation = Column(String, index=True)
    annotation_highest_impact = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class LowQualityClonalVariant(Base):

    __tablename__ = LOW_QUALITY_CLONAL_VARIANT_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
    overlaps_multiple_genes = Column(Boolean, default=False)

    annotation = Column(String, index=True)
    annotation_highest_impact = Column(String, index=True)
    annotation_impact = Column(String, index=True)
    gene_name = Column(String, index=True)
    gene_id = Column(String)
    biotype = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)
    cdna_pos_length = Column(String)
    cds_pos_length = Column(String)
    aa_pos_length = Column(String)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    def get_variant_id(self):
        return "{}:{}>{}".format(self.position, self.reference, self.alternate)


class VariantObservation(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = VARIANT_OBSERVATION_TABLE_NAME

    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    quality = Column(Float)
    filter = Column(String)
    """
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
    ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
    ##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
    """
    dp = Column(Integer)
    ac = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    vaf = Column(Float)
    strand_bias = Column(Integer)
    # fields replicated from Variant for performance reasons
    annotation = Column(String)
    annotation_impact = Column(String)
    biotype = Column(String)
    annotation_highest_impact = Column(String)
    gene_name = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)

    # fields replicated from sample for performance reasons
    date = Column(Date)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    ForeignKeyConstraint([sample], [SampleEna.run_accession])
    ForeignKeyConstraint([variant_id], [Variant.variant_id])

    __table_args__ = (
        Index("{}_i1".format(VARIANT_OBSERVATION_TABLE_NAME),
              "annotation_highest_impact", "position"),
        Index("{}_i2".format(VARIANT_OBSERVATION_TABLE_NAME), "sample"),
        Index("{}_i3".format(VARIANT_OBSERVATION_TABLE_NAME), "position"),
        Index("{}_i4".format(VARIANT_OBSERVATION_TABLE_NAME), "annotation_highest_impact"),
        Index("{}_i5".format(VARIANT_OBSERVATION_TABLE_NAME), "variant_id")
    )


class VariantObservationCovid19Portal(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME

    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    quality = Column(Float)
    filter = Column(String)
    """
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
    ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    ##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">
    ##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">
    """
    dp = Column(Integer)
    ac = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    vaf = Column(Float)
    strand_bias = Column(Integer)
    # fields replicated from Variant for performance reasons
    annotation = Column(String)
    annotation_impact = Column(String)
    biotype = Column(String)
    annotation_highest_impact = Column(String)
    gene_name = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)

    # fields replicated from sample for performance reasons
    date = Column(Date)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    ForeignKeyConstraint([sample], [SampleCovid19Portal.run_accession])
    ForeignKeyConstraint([variant_id], [VariantCovid19Portal.variant_id])

    __table_args__ = (
        Index("{}_i1".format(COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME),
              "annotation_highest_impact", "position"),
        Index("{}_i2".format(COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME), "sample"),
        Index("{}_i3".format(COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME), "position"),
        Index("{}_i4".format(COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME), "annotation_highest_impact"),
        Index("{}_i5".format(COVID19_PORTAL_VARIANT_OBSERVATION_TABLE_NAME), "variant_id")
    )


class SubclonalVariantObservation(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = SUBCLONAL_VARIANT_OBSERVATION_TABLE_NAME

    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    quality = Column(Float)
    filter = Column(String)
    dp = Column(Integer)
    ac = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    vaf = Column(Float)
    strand_bias = Column(Integer)

    # fields replicated from Variant for performance reasons
    annotation = Column(String)
    annotation_impact = Column(String)
    biotype = Column(String)
    annotation_highest_impact = Column(String, index=True)
    gene_name = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)

    # fields replicated from sample for performance reasons
    date = Column(Date)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    ForeignKeyConstraint([sample], [SampleEna.run_accession])
    ForeignKeyConstraint([variant_id], [SubclonalVariant.variant_id])

    __table_args__ = (
        Index("{}_i1".format(SUBCLONAL_VARIANT_OBSERVATION_TABLE_NAME), "position"),
        Index("{}_i2".format(SUBCLONAL_VARIANT_OBSERVATION_TABLE_NAME), "annotation_highest_impact", "vaf"),
        Index("{}_i3".format(SUBCLONAL_VARIANT_OBSERVATION_TABLE_NAME), "vaf"),
    )


class LowFrequencyVariantObservation(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = LOW_FREQUENCY_VARIANT_OBSERVATION_TABLE_NAME

    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    quality = Column(Float)
    filter = Column(String)
    dp = Column(Integer)
    ac = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    vaf = Column(Float)
    strand_bias = Column(Integer)

    # fields replicated from Variant for performance reasons
    annotation = Column(String)
    annotation_impact = Column(String)
    biotype = Column(String)
    annotation_highest_impact = Column(String, index=True)
    gene_name = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)

    # fields replicated from sample for performance reasons
    date = Column(Date)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    ForeignKeyConstraint([sample], [SampleEna.run_accession])
    ForeignKeyConstraint([variant_id], [LowFrequencyVariant.variant_id])

    __table_args__ = (
        Index("{}_i1".format(LOW_FREQUENCY_VARIANT_OBSERVATION_TABLE_NAME), "position"),
        Index("{}_i2".format(LOW_FREQUENCY_VARIANT_OBSERVATION_TABLE_NAME), "annotation_highest_impact", "vaf"),
        Index("{}_i3".format(LOW_FREQUENCY_VARIANT_OBSERVATION_TABLE_NAME), "vaf"),
    )


class LowQualityClonalVariantObservation(Base):
    """
    A variant observation in a particular sample. This contains all annotations of a specific observation of a variant.
    """
    __tablename__ = LOW_QUALITY_CLONAL_VARIANT_OBSERVATION_TABLE_NAME

    sample = Column(String, primary_key=True)
    variant_id = Column(String, primary_key=True)
    chromosome = Column(String)
    position = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    quality = Column(Float)
    filter = Column(String)
    dp = Column(Integer)
    ac = Column(Integer)
    dp4_ref_forward = Column(Integer)
    dp4_ref_reverse = Column(Integer)
    dp4_alt_forward = Column(Integer)
    dp4_alt_reverse = Column(Integer)
    vaf = Column(Float)
    strand_bias = Column(Integer)

    # fields replicated from Variant for performance reasons
    annotation = Column(String)
    annotation_impact = Column(String)
    biotype = Column(String)
    annotation_highest_impact = Column(String, index=True)
    gene_name = Column(String)
    hgvs_c = Column(String)
    hgvs_p = Column(String)

    # fields replicated from sample for performance reasons
    date = Column(Date)

    # derived annotations
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    length = Column(Integer)
    reference_amino_acid = Column(String)
    alternate_amino_acid = Column(String)
    position_amino_acid = Column(Integer)

    # ConsHMM conservation annotations
    cons_hmm_sars_cov_2 = Column(Float)
    cons_hmm_sarbecovirus = Column(Float)
    cons_hmm_vertebrate_cov = Column(Float)

    # Pfam protein domains
    pfam_name = Column(String)
    pfam_description = Column(String)

    ForeignKeyConstraint([sample], [SampleEna.run_accession])
    ForeignKeyConstraint([variant_id], [LowQualityClonalVariant.variant_id])

    __table_args__ = (
        Index("{}_i1".format(LOW_QUALITY_CLONAL_VARIANT_OBSERVATION_TABLE_NAME), "position"),
        Index("{}_i2".format(LOW_QUALITY_CLONAL_VARIANT_OBSERVATION_TABLE_NAME), "annotation_highest_impact", "vaf"),
        Index("{}_i3".format(LOW_QUALITY_CLONAL_VARIANT_OBSERVATION_TABLE_NAME), "vaf"),
    )


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
    DOWNLOADER = 3


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


class LastUpdate(Base):

    __tablename__ = LAST_UPDATE_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    update_time = Column(Date(), nullable=False)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__), nullable=False)


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


class PrecomputedVariantsPerSample(Base):

    __tablename__ = PRECOMPUTED_VARIANTS_PER_SAMPLE_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    count = Column(Integer)
    number_mutations = Column(Integer)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    gene_name = Column(String)


class PrecomputedSubstitutionsCounts(Base):

    __tablename__ = PRECOMPUTED_SUBSTITUTIONS_COUNTS_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    count = Column(Integer)
    reference = Column(String)
    alternate = Column(String)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    gene_name = Column(String)


class PrecomputedIndelLength(Base):

    __tablename__ = PRECOMPUTED_INDEL_LENGTH_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    count = Column(Integer)
    length = Column(Integer)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))
    gene_name = Column(String)


class PrecomputedAnnotation(Base):

    __tablename__ = PRECOMPUTED_ANNOTATION_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    count = Column(Integer)
    annotation = Column(String)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))
    gene_name = Column(String)


class PrecomputedOccurrence(Base):

    __tablename__ = PRECOMPUTED_OCCURRENCE_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    total = Column(Integer)
    count = Column(Integer)
    month = Column(String)
    frequency = Column(Float)
    frequency_by_month = Column(Float)
    variant_id = Column(String)
    hgvs_p = Column(String)
    gene_name = Column(String)
    domain = Column(String)
    annotation = Column(String)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))


class RegionType(enum.Enum):

    __constraint_name__ = REGION_TYPE_CONSTRAINT_NAME

    GENE = 1
    DOMAIN = 2
    CODING_REGION = 3


class PrecomputedSynonymousNonSynonymousCounts(Base):

    __tablename__ = PRECOMPUTED_NS_S_COUNTS_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    month = Column(Date)
    region_type = Column(Enum(RegionType, name=RegionType.__constraint_name__))
    region_name = Column(String)
    country = Column(String)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))
    ns = Column(Integer)
    s = Column(Integer)


class PrecomputedTableCounts(Base):

    __tablename__ = PRECOMPUTED_TABLE_COUNTS_TABLE_NAME
    VIRTUAL_TABLE_COUNTRY = "Country"
    FACTOR_SOURCE = "source"

    id = Column(Integer, primary_key=True, autoincrement=True)
    table = Column(String)
    factor = Column(String)
    value = Column(String)
    count = Column(Integer)


class PrecomputedVariantAbundanceHistogram(Base):

    __tablename__ = PRECOMPUTED_VARIANT_ABUNDANCE_HIST_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    position_bin = Column(Integer)
    count_unique_variants = Column(Integer)
    count_variant_observations = Column(Integer)
    bin_size = Column(Integer)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))


class PrecomputedVariantsPerLineage(Base):

    __tablename__ = PRECOMPUTED_VARIANTS_PER_LINEAGE_TABLE_NAME

    id = Column(Integer, primary_key=True, autoincrement=True)
    lineage = Column(String)
    variant_id = Column(String)
    country = Column(String)
    count_observations = Column(Integer)
    source = Column(Enum(DataSource, name=DataSource.__constraint_name__))


class Lineages(Base):
    """
    Annotate pangolin lineage identifiers with WHO designation, VOC/VUI, parent name
    """
    __tablename__ = LINEAGE_TABLE_NAME

    pango_lineage_id = Column(String, primary_key=True)
    # Constellation label used by scorpio for assignment
    constellation_id = Column(String)
    who_label = Column(String)
    # VOC/VUI/V information
    phe_label = Column(String)
    voc_date = Column(Date)
    vui_date = Column(Date)
    variant_of_concern = Column(Boolean, default=False)
    variant_under_investigation = Column(Boolean, default=False)
    parent_lineage_id = Column(String)
    tags = Column(String)

    variants = relationship(
        "LineageDefiningVariants",
        secondary=LINEAGE_SITES_JUNCTION_TABLE_NAME,
        back_populates='constellations',
    )


class VariantLevel(enum.Enum):
    __constraint_name__ = VARIANT_LEVEL_CONSTRAINT_NAME

    PROTEOMIC = 1
    GENOMIC = 2


class LineageDefiningVariants(Base):
    """
    Store lineage defining mutations as defined in the scorpio constellation files. Mutations can be in nucleotide
    (intergenic) or proteomic space. The columns variant_id, position, reference and alternate therefore have different
    meanings and formatting

    Variant_id: position_in_genome:[reference_bases]>[alternate_nuc] (genomic level)
                canonical_protein_name:[reference_aa]position_in_protein[alternate_aa] (proteomic level)
    position: position_in_genome (genomic level) or position_in_protein (proteomic level)
    reference: reference_bases (genomic level) or reference_aa (proteomic level)
    alternate: alternate_bases (genomic level) or alternate_aa (proteomic level)
    """
    __tablename__ = CONSTELLATION_SITES_TABLE_NAME

    variant_id = Column(String, primary_key=True)
    variant_type = Column(Enum(VariantType, name=VariantType.__constraint_name__))
    protein = Column(String)
    position = Column(Integer, index=True)
    reference = Column(String)
    alternate = Column(String)
    ambiguous_alternate = Column(Boolean, default=False)
    annotation = Column(String)
    hgvs = Column(String)
    variant_level = Column(Enum(VariantLevel, name=VariantLevel.__constraint_name__))

    constellations = relationship(
        "Lineages",
        secondary=LINEAGE_SITES_JUNCTION_TABLE_NAME,
        back_populates="variants")


class LineageVariant(Base):
    """
    Junction table that maps constellation sites to their respective lineage
    """
    __tablename__ = LINEAGE_SITES_JUNCTION_TABLE_NAME

    pango_lineage_id = Column(ForeignKey(Lineages.pango_lineage_id), primary_key=True)
    variant_id = Column(ForeignKey(LineageDefiningVariants.variant_id), primary_key=True)


class NewsType(enum.Enum):
    __constraint_name__ = NEWS_TYPE_CONSTRAINT_NAME

    RELEASE = 1
    WARNINGS = 2
    BUGS = 3


class NewsSection(Base):
    """
    Table to hold information for news section
    """
    __tablename__ = NEWS_SECTION_TABLE_NAME

    message_text = Column(String, primary_key=True)
    published_date = Column(DateTime(timezone=True), nullable=False, default=datetime.now())
    message_type = Column(Enum(NewsType, name=NewsType.__constraint_name__))

