
-- The date stored in the variant observations tables need to be patched with the collection date from the corresponding sample.

-- patch ENA dataset
update variant_observation_on set date=s.collection_date from sample_ena_on as s where s.run_accession = variant_observation_on.sample;
update subclonal_variant_observation_on set date=s.collection_date from sample_ena_on as s where s.run_accession = subclonal_variant_observation_on.sample;
update low_frequency_variant_observation_on set date=s.collection_date from sample_ena_on as s where s.run_accession = low_frequency_variant_observation_on.sample;
update lq_clonal_variant_observation_on set date=s.collection_date from sample_ena_on as s where s.run_accession = lq_clonal_variant_observation_on.sample;

-- patch COVID19 Data Portal dataset
update variant_observation_covid19portal_on set date=s.collection_date from sample_covid19_portal_on as s where s.run_accession = variant_observation_covid19portal_on.sample;
