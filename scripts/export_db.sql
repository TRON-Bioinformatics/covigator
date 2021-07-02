-- this is exporting the ENA related tables
\copy variant_v11 to program 'gzip > variant_v11.csv.gz' csv header;
\copy variant_observation_v11 to program 'gzip > variant_observation_v11.csv.gz' csv header;
\copy sample_ena_v11 to program 'gzip > sample_ena_v11.csv.gz' csv header;
\copy sample_v11 to program 'gzip > sample_v11.csv.gz' csv header;
\copy job_ena_v11(run_accession,status,created_at,queued_at,downloaded_at,analysed_at,cleaned_at,loaded_at,cooccurrence_at,failed_at,error_message,fastq_path,vcf_path) to program 'gzip > job_ena_v11.csv.gz' csv header;
\copy subclonal_variant_observation_v11 to program 'gzip > subclonal_variant_observation_v11.csv.gz' csv header;
\copy variant_cooccurrence_v11 to program 'gzip > variant_cooccurrence_v11.csv.gz' csv header;

-- this is exporting the GISAID related tables
\copy sample_v12 to program 'gzip > sample_v12.csv.gz' csv header;
-- in this table we need to specifiy the columns to avoid exporting the table with the sequences which is large
\copy sample_gisaid_v12(run_accession, date, host_tax_id, host, country_raw, region, country, country_alpha_2, country_alpha_3, continent, continent_alpha_2, site, site2) to program 'gzip > sample_gisaid_v12.csv.gz' csv header;
\copy job_gisaid_v12 to program 'gzip > job_gisaid_v12.csv.gz' csv header;
\copy variant_observation_v12 to program 'gzip > variant_observation_v12.csv.gz' csv header;
\copy variant_v12 to program 'gzip > variant_v12.csv.gz' csv header;


