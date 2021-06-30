

ALTER TABLE variant_v13 ADD variant_type variant_type_v13;
ALTER TABLE variant_observation_v13 ADD variant_type variant_type_v13;
ALTER TABLE subclonal_variant_observation_v13 ADD variant_type variant_type_v13;

update variant_observation_v13 set variant_type='SNV'
    where length(reference)=1 and length(alternate)=1;
update variant_observation_v13 set variant_type='DELETION'
    where length(reference)>1 and length(alternate)=1;
update variant_observation_v13 set variant_type='INSERTION'
    where length(reference)=1 and length(alternate)>1;

update subclonal_variant_observation_v13 set variant_type='SNV'
    where length(reference)=1 and length(alternate)=1;
update subclonal_variant_observation_v13 set variant_type='DELETION'
    where length(reference)>1 and length(alternate)=1;
update subclonal_variant_observation_v13 set variant_type='INSERTION'
    where length(reference)=1 and length(alternate)>1;
