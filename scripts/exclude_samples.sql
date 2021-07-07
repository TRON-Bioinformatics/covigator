
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

-- exclude samples with a ratio of Ns or ambiguous bases greater than 0.2 or a sequence length < 20 % of the genome
update job_gisaid_v14 set status='EXCLUDED' where run_accession in (
    select run_accession from (
        select run_accession,
            cast(sequence_length as float) / 29903 as coverage_ratio,
            cast(count_n_bases + count_ambiguous_bases as float) / sequence_length as n_ratio
        from sample_gisaid_v14)
    as counts where n_ratio > 0.2 or coverage_ratio < 0.2);
