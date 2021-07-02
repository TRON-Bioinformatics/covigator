

ALTER TABLE variant_observation_v13 ADD column annotation_highest_impact varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column annotation_highest_impact varchar;
ALTER TABLE variant_v13 ADD column annotation_highest_impact varchar;

update variant_v13 set annotation_highest_impact=regexp_replace(annotation, '&.*', '');
update variant_observation_v13 set annotation_highest_impact=regexp_replace(annotation, '&.*', '');
update subclonal_variant_observation_v13 set annotation_highest_impact=regexp_replace(annotation, '&.*', '');
