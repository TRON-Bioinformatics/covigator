

ALTER TABLE variant_v13 ADD column cons_hmm_sars_cov_2 float;
ALTER TABLE variant_v13 ADD column cons_hmm_sarbecovirus float;
ALTER TABLE variant_v13 ADD column cons_hmm_vertebrate_cov float;
ALTER TABLE variant_v13 ADD column pfam_name varchar;
ALTER TABLE variant_v13 ADD column pfam_description varchar;

ALTER TABLE variant_observation_v13 ADD column cons_hmm_sars_cov_2 float;
ALTER TABLE variant_observation_v13 ADD column cons_hmm_sarbecovirus float;
ALTER TABLE variant_observation_v13 ADD column cons_hmm_vertebrate_cov float;
ALTER TABLE variant_observation_v13 ADD column pfam_name varchar;
ALTER TABLE variant_observation_v13 ADD column pfam_description varchar;
ALTER TABLE variant_observation_v13 ADD column annotation_impact varchar;
ALTER TABLE variant_observation_v13 ADD column biotype varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column cons_hmm_vertebrate_cov float;
ALTER TABLE subclonal_variant_observation_v13 ADD column pfam_name varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column pfam_description varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column annotation_impact varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column biotype varchar;

update variant_observation_v13 set annotation_impact=v.annotation_impact, biotype=v.biotype
    from variant_observation_v13 as vo join variant_v13 as v on v.variant_id=vo.variant_id;
update subclonal_variant_observation_v13 set annotation_impact=v.annotation_impact, biotype=v.biotype
    from subclonal_variant_observation_v13 as vo join variant_v13 as v on v.variant_id=vo.variant_id;

-- TODO: update the conservation values from conservation table

