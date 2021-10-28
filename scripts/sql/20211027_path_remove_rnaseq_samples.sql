

delete from variant_observation_v19 where sample in (select run_accession from sample_ena_v19 where library_strategy='RNA-Seq');
delete from subclonal_variant_observation_v19 where sample in (select run_accession from sample_ena_v19 where library_strategy='RNA-Seq');
update job_ena_v19 set status='EXCLUDED' where run_accession in (select run_accession from sample_ena_v19 where library_strategy='RNA-Seq');
update sample_ena_v19 set finished=false where run_accession in (select run_accession from job_ena_v19 where status!='FINISHED');