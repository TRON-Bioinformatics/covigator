#!/bin/bash


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

job_ena="\\copy job_ena$version(run_accession,status,created_at,queued_at,downloaded_at,analysed_at,cleaned_at,loaded_at,cooccurrence_at,failed_at,error_message,fastq_path,vcf_path,qc_path,horizontal_coverage_path,vertical_coverage_path,num_reads,covered_bases,coverage,mean_depth,mean_base_quality,mean_mapping_quality) from program 'gzip -dc $input_folder/job_ena.csv.gz' csv header;"
job_gisaid=`get_import_command "job_gisaid"`
log=`get_import_command "log"`

precomputed_annotation_delete=`get_delete_command "precomputed_annotation"`
precomputed_indel_length_delete=`get_delete_command "precomputed_indel_length"`
precomputed_substitutions_counts_delete=`get_delete_command "precomputed_substitutions_counts"`
precomputed_top_occurrence_delete=`get_delete_command "precomputed_top_occurrence"`
precomputed_variants_per_sample_delete=`get_delete_command "precomputed_variants_per_sample"`
precomputed_table_counts_delete=`get_delete_command "precomputed_table_counts"`
precomputed_variant_abundance_histogram_delete=`get_delete_command "precomputed_variant_abundance_histogram"`
precomputed_ns_s_counts_delete=`get_delete_command "precomputed_ns_s_counts"`

precomputed_annotation=`get_import_command "precomputed_annotation"`
precomputed_indel_length=`get_import_command "precomputed_indel_length"`
precomputed_substitutions_counts=`get_import_command "precomputed_substitutions_counts"`
precomputed_top_occurrence=`get_import_command "precomputed_top_occurrence"`
precomputed_variants_per_sample=`get_import_command "precomputed_variants_per_sample"`
precomputed_table_counts=`get_import_command "precomputed_table_counts"`
precomputed_variant_abundance_histogram=`get_import_command "precomputed_variant_abundance_histogram"`
precomputed_ns_s_counts=`get_import_command "precomputed_ns_s_counts"`

sample_ena=`get_import_command "sample_ena"`
# NOTE: add column finished back in
sample_gisaid="\\copy sample_gisaid$version(run_accession,date,host_tax_id,host,country_raw,region,country,country_alpha_2,country_alpha_3,continent,continent_alpha_2,site,site2,sequence_length,count_n_bases,count_ambiguous_bases,count_snvs,count_insertions,count_deletions,finished) from program 'gzip -dc $input_folder/sample_gisaid.csv.gz' csv header;"
sample=`get_import_command "sample"`
subclonal_variant_observation=`get_import_command "subclonal_variant_observation"`
variant_cooccurrence=`get_import_command "variant_cooccurrence"`
variant_observation=`get_import_command "variant_observation"`
variant=`get_import_command "variant"`

psql $pg_uri -c "$sample_ena"
psql $pg_uri -c "$job_ena"
psql $pg_uri -c "$sample_gisaid"
psql $pg_uri -c "$job_gisaid"
psql $pg_uri -c "$sample"
psql $pg_uri -c "$variant"
psql $pg_uri -c "$subclonal_variant_observation"
psql $pg_uri -c "$variant_cooccurrence"
psql $pg_uri -c "$variant_observation"
psql $pg_uri -c "$log"

psql $pg_uri -c "$precomputed_annotation_delete"
psql $pg_uri -c "$precomputed_annotation"
psql $pg_uri -c "$precomputed_indel_length_delete"
psql $pg_uri -c "$precomputed_indel_length"
psql $pg_uri -c "$precomputed_substitutions_counts_delete"
psql $pg_uri -c "$precomputed_substitutions_counts"
psql $pg_uri -c "$precomputed_top_occurrence_delete"
psql $pg_uri -c "$precomputed_top_occurrence"
psql $pg_uri -c "$precomputed_variants_per_sample_delete"
psql $pg_uri -c "$precomputed_variants_per_sample"
psql $pg_uri -c "$precomputed_table_counts_delete"
psql $pg_uri -c "$precomputed_table_counts"
psql $pg_uri -c "$precomputed_variant_abundance_histogram_delete"
psql $pg_uri -c "$precomputed_variant_abundance_histogram"
#psql $pg_uri -c "$precomputed_ns_s_counts_delete"
#psql $pg_uri -c "$precomputed_ns_s_counts"
