#!/bin/bash


# USAGE: import_db.sh covigator_config.txt /your/data/folder

source $1
input_folder=$2

version=$COVIGATOR_TABLE_VERSION
export PGPASSWORD=$COVIGATOR_DB_PASSWORD
pg_uri=postgresql://$COVIGATOR_DB_USER@$COVIGATOR_DB_HOST:$COVIGATOR_DB_PORT/$COVIGATOR_DB_NAME


get_import_command() {
  echo "\\copy $1$version from program 'gzip -dc $input_folder/$1.csv.gz' csv header;"
}

get_delete_command() {
  echo "delete from $1$version;"
}

load_table() {
  psql $pg_uri -c "`get_delete_command $1`"
  psql $pg_uri -c "`get_import_command $1`"
}

# references
load_table "conservation"
load_table "gene"
load_table "domain"

# logs
load_table "log"
load_table "last_update"

# precomputations
load_table "precomputed_annotation"
load_table "precomputed_indel_length"
load_table "precomputed_ns_s_counts"
load_table "precomputed_substitutions_counts"
load_table "precomputed_table_counts"
load_table "precomputed_top_occurrence"
load_table "precomputed_variant_abundance_histogram"
load_table "precomputed_variants_per_lineage"
load_table "precomputed_variants_per_sample"

# ENA
sample_ena_fields="run_accession, finished, sample_accession, scientific_name, study_accession, experiment_accession, \
first_created, collection_date, instrument_platform, instrument_model, library_name, nominal_length, \
library_layout, library_strategy, library_source, library_selection, read_count, base_count, \
sample_collection, sequencing_method, center_name, fastq_ftp, fastq_md5, num_fastqs, host_tax_id, \
host_sex, host_body_site, host_gravidity, host_phenotype, host_genotype, lat, lon, country_raw, country, \
country_alpha_2, country_alpha_3, continent, continent_alpha_2, count_snvs, count_insertions, count_deletions, \
count_subclonal_snvs, count_subclonal_insertions, count_subclonal_deletions, count_low_frequency_snvs, \
count_low_frequency_insertions, count_low_frequency_deletions, status, created_at, queued_at, downloaded_at, \
analysed_at, cleaned_at, loaded_at, cooccurrence_at, failed_at, error_message, sample_folder, fastq_path, \
lofreq_vcf_path, ivar_vcf_path, gatk_vcf_path, bcftools_vcf_path, lofreq_pangolin_path, ivar_pangolin_path, \
gatk_pangolin_path, bcftools_pangolin_path, fastp_path, deduplication_metrics_path, horizontal_coverage_path, \
vertical_coverage_path, num_reads, covered_bases, coverage, mean_depth, mean_base_quality, \
mean_mapping_quality, pangolin_lineage, pangolin_conflict, pangolin_ambiguity_score, pangolin_scorpio_call, \
pangolin_scorpio_support, pangolin_scorpio_conflict, pangolin_version, pangolin_pangolin_version, \
pangolin_pango_version, pangolin_status, pangolin_note, percent_duplication, \
unpaired_reads_examined, read_pairs_examined, secondary_or_supplementary_reads, unmapped_reads, \
unpaired_read_duplicates, read_pair_duplicates, read_pair_optical_duplicates, covigator_accessor_version, \
covigator_processor_version"
psql $pg_uri -c "\\copy sample_ena$version($sample_ena_fields) from program 'gzip -dc $input_folder/sample_ena.csv.gz' csv header;"
load_table "variant"
load_table "variant_cooccurrence"
load_table "variant_observation"
load_table "subclonal_variant"
load_table "subclonal_variant_observation"
load_table "low_frequency_variant_observation"

# GISAID
sample_gisaid_fields="run_accession, finished, collection_date, host_tax_id, host, \
 country_raw, region, country, country_alpha_2, country_alpha_3, continent, continent_alpha_2, \
 site, site2, sequence_length, count_n_bases, count_ambiguous_bases, count_snvs, \
 count_insertions, count_deletions, status, created_at, queued_at, analysed_at, loaded_at, failed_at, \
 error_message, sample_folder, vcf_path, fasta_path, pangolin_path, pangolin_lineage, \
 pangolin_conflict, pangolin_ambiguity_score, pangolin_scorpio_call, pangolin_scorpio_support, \
 pangolin_scorpio_conflict, pangolin_version, pangolin_pangolin_version, \
 pangolin_pango_version, pangolin_status, pangolin_note, covigator_accessor_version, covigator_processor_version"
psql $pg_uri -c "\\copy sample_gisaid$version($sample_gisaid_fields) from program 'gzip -dc $input_folder/sample_gisaid.csv.gz' csv header;"
load_table "gisaid_variant"
load_table "gisaid_variant_observation"
load_table "gisaid_variant_cooccurrence"
