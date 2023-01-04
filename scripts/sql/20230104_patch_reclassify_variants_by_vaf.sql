
-- move variant observations with VAF >= 0.5 and < 0.8 to the low quality clonal tables
insert into lq_clonal_variant_on (select * from subclonal_variant_on where variant_id in (select variant_id from subclonal_variant_observation_on where vaf >= 0.5));
insert into lq_clonal_variant_observation_on (select * from subclonal_variant_observation_on where vaf >= 0.5);
delete from subclonal_variant_observation_on where vaf >= 0.5;
--delete from subclonal_variant_on where variant_id not in (select variant_id from subclonal_variant_observation_on);


-- move low frequency variant observations with VAF >= 0.02 to the subclonal tables
insert into subclonal_variant_on (select * from low_frequency_variant_on where variant_id in (
    select variant_id from low_frequency_variant_observation_on where vaf >= 0.02)) on conflict do nothing;
insert into subclonal_variant_observation_on (select * from low_frequency_variant_observation_on where vaf >= 0.02);
delete from low_frequency_variant_observation_on where vaf >= 0.02;