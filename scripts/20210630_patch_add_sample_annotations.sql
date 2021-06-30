

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
