-- size of whole database in disk
SELECT pg_size_pretty (pg_database_size ('covigator'));

-- size of table
select pg_size_pretty(pg_relation_size('variant_observation_v12'));

-- size of indexes
SELECT pg_size_pretty(pg_indexes_size('variant_observation_v12'));

-- percentage of finished jobs
select (select count(*) from job_gisaid_v12 where status='FINISHED' or status='FAILED_PROCESSING')::float/
(select count(*) from job_gisaid_v12);