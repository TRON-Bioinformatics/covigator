
CREATE INDEX ON sample_gisaid_v15(date, country, finished);
CREATE INDEX ON sample_gisaid_v15(date, country);
CREATE INDEX ON sample_ena_v15(collection_date, country, finished);
CREATE INDEX ON sample_ena_v15(collection_date, country);

CREATE INDEX ON variant_observation_v14(variant_id, date_trunc('month', date::timestamp));
CREATE INDEX ON variant_observation_v14(variant_id, hgvs_p, gene_name, annotation_highest_impact);