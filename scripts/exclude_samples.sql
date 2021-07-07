
-- exclude samples from a too early date
update job_gisaid_v12 set status='EXCLUDED' where run_accession in (
    select run_accession from sample_gisaid_v12 where date < '2019-12-01');

-- exclude samples with more than 76 SNVs
update job_gisaid_v12 set status='EXCLUDED' where run_accession in (
    select run_accession from (select count(*) as count, sample from variant_observation_v12
    where length(reference) = 1 and length(alternate) = 1 group by sample) as counts where count > 76);

-- exclude samples with more than 10 deletions
update job_gisaid_v12 set status='EXCLUDED' where run_accession in (
    select run_accession from (select count(*) as count, sample from variant_observation_v12
    where length(reference) > 1 and length(alternate) = 1 group by sample) as counts where count > 10);

-- exclude samples with more than 10 insertions
update job_gisaid_v12 set status='EXCLUDED' where run_accession in (
    select run_accession from (select count(*) as count, sample from variant_observation_v12
    where length(reference) = 1 and length(alternate) > 1 group by sample) as counts where count > 10);
