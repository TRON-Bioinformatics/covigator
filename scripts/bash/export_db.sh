#!/bin/bash


# USAGE: import_db.sh covigator_config.txt /your/data/folder

source $1
output=$2

version=$COVIGATOR_TABLE_VERSION
pg_uri=postgresql://$COVIGATOR_DB_USER:$COVIGATOR_DB_PASSWORD@$COVIGATOR_DB_HOST:$COVIGATOR_DB_PORT/$COVIGATOR_DB_NAME

mkdir -p $output

get_export_command() {
  echo "\\copy $1$version to program 'gzip > $output/$1.csv.gz' csv header;"
}

# references
psql $pg_uri -c "$(get_export_command "conservation")"
psql $pg_uri -c "$(get_export_command "gene")"
psql $pg_uri -c "$(get_export_command "domain")"

# logs
psql $pg_uri -c "$(get_export_command "log")"
psql $pg_uri -c "$(get_export_command "last_update")"

# precomputations
psql $pg_uri -c "$(get_export_command "precomputed_annotation")"
psql $pg_uri -c "$(get_export_command "precomputed_indel_length")"
psql $pg_uri -c "$(get_export_command "precomputed_ns_s_counts")"
psql $pg_uri -c "$(get_export_command "precomputed_substitutions_counts")"
psql $pg_uri -c "$(get_export_command "precomputed_table_counts")"
psql $pg_uri -c "$(get_export_command "precomputed_top_occurrence")"
psql $pg_uri -c "$(get_export_command "precomputed_variant_abundance_histogram")"
psql $pg_uri -c "$(get_export_command "precomputed_variants_per_lineage")"
psql $pg_uri -c "$(get_export_command "precomputed_variants_per_sample")"

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
pangolin_scorpio_version, pangolin_constellation_version, pangolin_qc_status, pangolin_qc_notes, pangolin_note, percent_duplication, \
unpaired_reads_examined, read_pairs_examined, secondary_or_supplementary_reads, unmapped_reads, \
unpaired_read_duplicates, read_pair_duplicates, read_pair_optical_duplicates, covigator_accessor_version, \
covigator_processor_version, intrahost_filter, potential_coinfection"
psql $pg_uri -c "\\copy sample_ena$version($sample_ena_fields) to program 'gzip > $output/sample_ena.csv.gz' csv header;"
psql $pg_uri -c "$(get_export_command "variant")"
psql $pg_uri -c "$(get_export_command "variant_cooccurrence")"
psql $pg_uri -c "$(get_export_command "variant_observation")"
psql $pg_uri -c "$(get_export_command "subclonal_variant")"
psql $pg_uri -c "$(get_export_command "subclonal_variant_observation")"
psql $pg_uri -c "$(get_export_command "low_frequency_variant")"
psql $pg_uri -c "$(get_export_command "low_frequency_variant_observation")"
psql $pg_uri -c "$(get_export_command "lq_clonal_variant")"
psql $pg_uri -c "$(get_export_command "lq_clonal_variant_observation")"

# Covid19portal
psql $pg_uri -c "$(get_export_command "sample_covid19_portal")"
psql $pg_uri -c "$(get_export_command "variant_covid19portal")"
psql $pg_uri -c "$(get_export_command "variant_observation_covid19portal")"
