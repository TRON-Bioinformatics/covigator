-- variants is the only table that can have overlapping entries for different samples
-- it requires an intermediate table to manage conflicts
\copy variant_v13 from program 'gzip -dc variant_v11.csv.gz' csv header;
CREATE TABLE tmp_variant_v13 (LIKE variant_v13 INCLUDING DEFAULTS);
\copy tmp_variant_v13 from program 'gzip -dc variant_v12.csv.gz' csv header;
INSERT INTO variant_v13 SELECT * FROM tmp_variant_v13 ON CONFLICT DO NOTHING;
drop table tmp_variant_v13;

\copy sample_ena_v13 from program 'gzip -dc sample_ena_v11.csv.gz' csv header;
\copy sample_gisaid_v13(run_accession, date, host_tax_id, host, country_raw, region, country, country_alpha_2, country_alpha_3, continent, continent_alpha_2, site, site2) from program 'gzip -dc sample_gisaid_v12.csv.gz' csv header;

\copy sample_v13 from program 'gzip -dc sample_v11.csv.gz' csv header;
\copy sample_v13 from program 'gzip -dc sample_v12.csv.gz' csv header;

\copy variant_observation_v13 from program 'gzip -dc variant_observation_v11.csv.gz' csv header;
\copy variant_observation_v13 from program 'gzip -dc variant_observation_v12.csv.gz' csv header;
\copy subclonal_variant_observation_v13 from program 'gzip -dc subclonal_variant_observation_v11.csv.gz' csv header;
\copy variant_cooccurrence_v11 from program 'gzip -dc variant_cooccurrence_v11.csv.gz' csv header;

\copy job_ena_v13(run_accession,status,created_at,queued_at,downloaded_at,analysed_at,cleaned_at,loaded_at,cooccurrence_at,failed_at,error_message,fastq_path,vcf_path) from program 'gzip -dc job_ena_v11.csv.gz' csv header;
\copy job_gisaid_v12 from program 'gzip -dc job_gisaid_v12.csv.gz' csv header;
