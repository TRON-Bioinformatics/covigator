

ALTER TABLE variant_observation_v13 ADD COLUMN annotation varchar;
ALTER TABLE variant_observation_v13 ADD COLUMN gene_name varchar;
ALTER TABLE variant_observation_v13 ADD COLUMN hgvs_p varchar;
ALTER TABLE variant_observation_v13 ADD COLUMN hgvs_c varchar;

ALTER TABLE subclonal_variant_observation_v13 ADD COLUMN annotation varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD COLUMN gene_name varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD COLUMN hgvs_p varchar;
ALTER TABLE subclonal_variant_observation_v13 ADD COLUMN hgvs_c varchar;

update variant_observation_v13 set annotation=v.annotation, gene_name=v.gene_name,
    hgvs_p=v.hgvs_p, hgvs_c=v.hgvs_c from variant_v13 as v where variant_observation_v13.variant_id=v.variant_id;
update subclonal_variant_observation_v13 set annotation=v.annotation, gene_name=v.gene_name,
    hgvs_p=v.hgvs_p, hgvs_c=v.hgvs_c from variant_v13 as v where subclonal_variant_observation_v13.variant_id=v.variant_id;


CREATE INDEX ON variant_observation_v13(variant_id, hgvs_p, gene_name, annotation);
