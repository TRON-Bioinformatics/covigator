

ALTER TABLE sample_gisaid_v13 ADD COLUMN finished boolean default False;
ALTER TABLE sample_ena_v13 ADD COLUMN finished boolean default False;


update sample_gisaid_v13 set finished=True where run_accession in (select run_accession from job_gisaid_v13 where status='FINISHED');
update sample_ena_v13 set finished=True where run_accession in (select run_accession from job_ena_v13 where status='FINISHED');

CREATE INDEX ON sample_gisaid_v13(date, country, finished);
CREATE INDEX ON sample_ena_v13(first_created, country, finished);
CREATE INDEX ON sample_ena_v13(first_created, country);