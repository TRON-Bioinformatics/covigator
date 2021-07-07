

ALTER TABLE variant_observation_v13 ADD column length integer;
ALTER TABLE subclonal_variant_observation_v13 ADD column length integer;
ALTER TABLE variant_v13 ADD column length integer;

update variant_v13 set length=length(alternate)-length(reference);
update variant_observation_v13 set length=length(alternate)-length(reference);
update subclonal_variant_observation_v13 set length=length(alternate)-length(reference);

ALTER TABLE variant_observation_v13 ADD column reference_amino_acid varchar;
ALTER TABLE variant_observation_v13 ADD column alternate_amino_acid varchar;
ALTER TABLE variant_observation_v13 ADD column position_amino_acid integer;
ALTER TABLE variant_v13 ADD column reference_amino_acid varchar;
ALTER TABLE variant_v13 ADD column alternate_amino_acid varchar;
ALTER TABLE variant_v13 ADD column position_amino_acid integer;
ALTER TABLE subclonal_variant_observation_v13 ADD column reference_amino_acid varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column alternate_amino_acid varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD column position_amino_acid integer;


