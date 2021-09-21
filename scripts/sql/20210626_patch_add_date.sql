

ALTER TABLE variant_observation_v13 ADD COLUMN date date;
ALTER TABLE subclonal_variant_observation_v13 ADD COLUMN date date;

-- this index is for the update to work better
CREATE INDEX ON variant_observation_v13(sample, source);

update variant_observation_v13 set date=s.first_created from sample_ena_v13 as s
    where variant_observation_v13.sample=s.run_accession and variant_observation_v13.source='ENA';
update variant_observation_v13 set date=s.date from sample_gisaid_v13 as s
    where variant_observation_v13.sample=s.run_accession and variant_observation_v13.source='GISAID';
update subclonal_variant_observation_v13 set date=s.first_created from sample_ena_v13 as s
    where subclonal_variant_observation_v13.sample=s.run_accession and subclonal_variant_observation_v13.source='ENA';

CREATE INDEX ON variant_observation_v13(variant_id, date_trunc('month', date::timestamp));
