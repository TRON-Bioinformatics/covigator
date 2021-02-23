SELECT p1.relname as table,
    p2.n_tup_ins - p2.n_tup_del as rows,
    pg_size_pretty(pg_relation_size(p1.relid)) as table_size,
    pg_size_pretty (pg_indexes_size(p1.relid)) as indexes_size
from pg_catalog.pg_statio_user_tables as p1,
    pg_stat_all_tables as p2 where p1.relid=p2.relid;