

ALTER TABLE sample_ena_v13 ADD column count_snvs integer;
ALTER TABLE sample_ena_v13 ADD column count_insertions integer;
ALTER TABLE sample_ena_v13 ADD column count_deletions integer;
ALTER TABLE sample_ena_v13 ADD column count_subclonal_snvs integer;
ALTER TABLE sample_ena_v13 ADD column count_subclonal_insertions integer;
ALTER TABLE sample_ena_v13 ADD column count_subclonal_deletions integer;
ALTER TABLE sample_gisaid_v13 ADD column count_snvs integer;
ALTER TABLE sample_gisaid_v13 ADD column count_insertions integer;
ALTER TABLE sample_gisaid_v13 ADD column count_deletions integer;
ALTER TABLE sample_gisaid_v13 ADD column sequence_length integer;
ALTER TABLE sample_gisaid_v13 ADD column count_n_bases integer;
ALTER TABLE sample_gisaid_v13 ADD column count_ambiguous_bases integer;


update sample_ena_v13 set count_snvs = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_ena_v13.run_accession and vo.variant_type='SNV');
update sample_ena_v13 set count_insertions = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_ena_v13.run_accession and vo.variant_type='INSERTION');
update sample_ena_v13 set count_deletions = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_ena_v13.run_accession and vo.variant_type='DELETION');
update sample_gisaid_v13 set count_snvs = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_gisaid_v13.run_accession and vo.variant_type='SNV');
update sample_gisaid_v13 set count_insertions = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_gisaid_v13.run_accession and vo.variant_type='INSERTION');
update sample_gisaid_v13 set count_deletions = (select count(*) from variant_observation_v13 as vo
    where vo.sample=sample_gisaid_v13.run_accession and vo.variant_type='DELETION');