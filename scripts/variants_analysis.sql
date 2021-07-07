-- count variants
select count(*) from variant_observation_v12 where length(reference) = 1 and length(alternate) = 1;
select count(*) from variant_observation_v12 where length(reference) > 1 and length(alternate) = 1;
select count(*) from variant_observation_v12 where length(reference) = 1 and length(alternate) > 1;

-- average, STD, min and max variants per sample
select avg(count), stddev(count), min(count), max(count) from (
    select count(*) as count, sample from variant_observation_v12 where length(reference) = 1 and length(alternate) = 1 group by sample
    ) as counts;

select avg(count), stddev(count), min(count), max(count) from (
    select count(*) as count, sample from variant_observation_v12 where length(reference) > 1 group by sample
    ) as counts;

select avg(count), stddev(count), min(count), max(count) from (
    select count(*) as count, sample from variant_observation_v12 where length(alternate) > 1 group by sample
    ) as counts;

-- 25, 50 and 75 percentiles of the number of variants per sample
select percentile_disc(0.25) within group (order by count),
       percentile_disc(0.5) within group (order by count),
       percentile_disc(0.75) within group (order by count) from (
    select count(*) as count, sample from variant_observation_v12
    where length(reference) = 1 and length(alternate) = 1 group by sample) as counts;

select percentile_disc(0.25) within group (order by count),
       percentile_disc(0.5) within group (order by count),
       percentile_disc(0.75) within group (order by count) from (
    select count(*) as count, sample from variant_observation_v12
    where length(reference) = 1 and length(alternate) > 1 group by sample) as counts;

select percentile_disc(0.25) within group (order by count),
       percentile_disc(0.5) within group (order by count),
       percentile_disc(0.75) within group (order by count) from (
    select count(*) as count, sample from variant_observation_v12
    where length(reference) > 1 and length(alternate) = 1 group by sample) as counts;

-- count samples with more than 92 variants
select count(*) from (
    select count(*) as count, sample from variant_observation_v12
    where variant_type='SNV' group by sample) as counts where count > 76;
select count(*) from (
    select count(*) as count, sample from variant_observation_v12
    where variant_type='DELETION' group by sample) as counts where count > 10;

-- indel length distribution
select count(*), length(reference) from variant_observation_v12 where length(reference) > 1 group by length(reference);
select count(*), length(alternate) from variant_observation_v12 where length(alternate) > 1 group by length(alternate);
\copy (select count(*), length(reference) from variant_observation_v12
    where length(reference) > 1 group by length(reference)) to 'length_distribution_deletion' csv header;
\copy (select count(*), length(alternate) from variant_observation_v12
    where length(alternate) > 1 group by length(alternate)) to 'length_distribution_insertion' csv header;
\copy (select count(*), length(reference) from variant_v12
    where length(reference) > 1 group by length(reference)) to 'length_distribution_deletion_unique.csv' csv header;
\copy (select count(*), length(alternate) from variant_v12
    where length(alternate) > 1 group by length(alternate)) to 'length_distribution_insertion_unique.csv' csv header;


-- compute dN/dS
select count(*), gene_name from variant_observation_v12 as vo join variant_v12 as v on vo.variant_id=v.variant_id
    where v.annotation = 'missense_variant' group by gene_name;
select count(*), gene_name from variant_observation_v12 as vo join variant_v12 as v on vo.variant_id=v.variant_id
    where v.annotation = 'synonymous_variant' group by gene_name;
\copy (select count(*) as dN, gene_name from variant_observation_v12 as vo
    join variant_v12 as v on vo.variant_id=v.variant_id
    where v.annotation = 'missense_variant' group by gene_name) to 'dn.csv' header csv;
\copy (select count(*) as dS, gene_name from variant_observation_v12 as vo
    join variant_v12 as v on vo.variant_id=v.variant_id
    where v.annotation = 'synonymous_variant' group by gene_name) to 'dn.csv' header csv;

-- distribution of aminoacid substitutions
select count(*), reference, alternate from variant_observation_v12
    where length(reference) = 1 and length(alternate) = 1 group by reference, alternate;

-- sample selection to test alignment settings
select * from (select count(*) as count, sample from variant_observation_v12
    where length(reference) > 8 and length(alternate) = 1 group by sample) as counts where count > 3 limit 10;
 count |              sample
-------+-----------------------------------
     6 | hCoV-19/Afghanistan/IMB07968/2020
     5 | hCoV-19/Albania/210901064/2021
     4 | hCoV-19/Albania/210901066/2021
     7 | hCoV-19/Albania/210901072/2021
     4 | hCoV-19/Angola/KRISP-K009681/2021
     5 | hCoV-19/Angola/KRISP-K009684/2021
     6 | hCoV-19/Angola/KRISP-K009688/2021
     4 | hCoV-19/Angola/KRISP-K009740/2021
     4 | hCoV-19/Angola/KRISP-K009769/2021
     4 | hCoV-19/Angola/KRISP-K009770/2021

select * from (select count(*) as count, sample from variant_observation_v12
    where length(alternate) > 8 and length(reference) = 1 group by sample) as counts where count > 3 limit 11;
 count |                 sample
-------+----------------------------------------
     5 | hCoV-19/Egypt/MASRI-C4-004/2020
     4 | hCoV-19/Egypt/MASRI-C5-010/2020
   195 | hCoV-19/Georgia/Tb-72720/2020
     5 | hCoV-19/Germany/NI-IOV-60240286/2021
     4 | hCoV-19/Germany/NI-IOV-60241720/2021
     4 | hCoV-19/Germany/NI-IOV-7786709585/2021
     5 | hCoV-19/Germany/NI-IOV-7790444785/2021
     4 | hCoV-19/Germany/NI-IOV-MHH_30/2020
     4 | hCoV-19/Germany/NI-RKI-I-014788/2020
     5 | hCoV-19/Germany/NW-RKI-I-014263/2021
     4 | hCoV-19/Germany/NW-RKI-I-019153/2021


-- FP calls from bad samples
select count(*) from variant_observation_v14 where run_accession in (
    select run_accession from sample_ena_v14 where mean_mapping_quality < 10);
select count(*) from variant_observation_v14 where sample in (
    select run_accession from job_ena_v14 where coverage < 20.0);

