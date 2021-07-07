-- size of whole database in disk
SELECT pg_size_pretty (pg_database_size ('covigator'));

-- size of all tables and indices and number of rows
SELECT p1.relname as table,
    p2.n_tup_ins - p2.n_tup_del as rows,
    pg_size_pretty(pg_relation_size(p1.relid)) as table_size,
    pg_size_pretty (pg_indexes_size(p1.relid)) as indexes_size
from pg_catalog.pg_statio_user_tables as p1,
    pg_stat_all_tables as p2 where p1.relid=p2.relid
order by pg_relation_size(p1.relid);

-- size of table
select pg_size_pretty(pg_relation_size('variant_observation_v12'));

-- size of indexes
SELECT pg_size_pretty(pg_indexes_size('variant_observation_v12'));

-- percentage of finished jobs
select (select count(*) from job_gisaid_v12 where status='FINISHED' or status='FAILED_PROCESSING')::float/
(select count(*) from job_gisaid_v12);