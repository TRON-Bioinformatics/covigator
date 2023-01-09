
-- apply filters on variants
insert into low_frequency_variant_on
    (select * from subclonal_variant_on where variant_id in
        (select variant_id from subclonal_variant_observation_on where ac < 10 or dp < 100 or length > 10 or length < -10)) on conflict do nothing;
insert into low_frequency_variant_observation_on
    (select * from subclonal_variant_observation_on where ac < 10 or dp < 100 or length > 10 or length < -10);
delete from subclonal_variant_observation_on where ac < 10 or dp < 100 or length > 10 or length < -10;
delete from subclonal_variant_on where variant_id not in (select distinct(variant_id) from subclonal_variant_observation_on);


-- add flag for low quality samples for intrahost calling to samples table
ALTER TABLE sample_ena_on ADD COLUMN intrahost_filter boolean default False;
-- flags low quality samples for intrahost calling
update sample_ena_on set intrahost_filter=true where read_count < 50000;
update sample_ena_on set intrahost_filter=true where covered_bases < 29000;
-- calculate the 0.999 percentile of intrahost mutations per sample
--select percentile_disc(0.999) within group (order by t.count_variants) from
--    (select s.c + l.c as count_variants, s.sample from
--        (select count(*) as c, sample from subclonal_variant_observation_on group by sample) as s,
--        (select count(*) as c, sample from low_frequency_variant_observation_on group by sample) as l
--        where s.sample = l.sample) as t;
-- this filters out samples with a number of mutations with VAF < 80 % over the 0.999 percentile
update sample_ena_on  set intrahost_filter=true from
    (select sum(s.c + l.c) as k, se.run_accession from sample_ena_on as se,
        (select sample, count(*) as c from subclonal_variant_observation_on group by sample) as s,
        (select sample, count(*) as c from low_frequency_variant_observation_on group by sample) as l
    where se.run_accession=s.sample and se.run_accession=l.sample group by se.run_accession) as t
where t.run_accession=sample_ena_on.run_accession and t.k >= 5147;


-- add flag for potential co-infection to samples table
ALTER TABLE sample_ena_on ADD COLUMN potential_coinfection boolean default False;
-- calculate the outlier threshold for mutations with VAF >= 40 % and VAF < 80 %
--select percentile_disc(0.5) + (2 * (percentile_disc(0.75) - percentile_disc(0.25))) within group (order by t.count_variants) from
--    (select s.c + l.c as count_variants, s.sample from
--        (select count(*) as c, sample from subclonal_variant_observation_on where vaf >= 0.4 group by sample) as s) as t;
-- flags potential co-infection to samples
update sample_ena_on  set potential_coinfection=true from
    (select sum(s.c) as k, se.run_accession from sample_ena_on as se,
        (select sample, count(*) as c from subclonal_variant_observation_on where vaf >= 0.4 group by sample) as s
    where se.run_accession=s.sample group by se.run_accession) as t
where t.run_accession=sample_ena_on.run_accession and t.k > 3;


insert into low_frequency_variant_on
    (select * from subclonal_variant_on where variant_id in
        (select variant_id from subclonal_variant_observation_on where sample in
            (select run_accession from sample_ena_on where intrahost_filter=true or potential_coinfection=true)
        )
    ) on conflict do nothing;
insert into low_frequency_variant_observation_on
    (select * from subclonal_variant_observation_on where sample in
        (select run_accession from sample_ena_on where intrahost_filter=true or potential_coinfection=true));
delete from subclonal_variant_observation_on where sample in
        (select run_accession from sample_ena_on where intrahost_filter=true or potential_coinfection=true);
delete from subclonal_variant_on where variant_id not in
    (select distinct(variant_id) from subclonal_variant_observation_on);
