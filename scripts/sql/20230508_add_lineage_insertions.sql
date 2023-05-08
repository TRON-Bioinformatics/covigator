INSERT INTO lineage_defining_variant_on(variant_id, variant_type, protein, position, reference, alternate, ambiguous_alternate, annotation, hgvs, variant_level) VALUES ('S:RD214REPED', 'INSERTION', 'S', '214', 'RD', 'REPED', FALSE, 'missense', 'p.R214_D215insEPE', 'aa');
INSERT INTO lineage_variant_on(pango_lineage_id, variant_id) VALUES('BA.1', 'S:RD214REPED');

INSERT INTO lineage_defining_variant_on(variant_id, variant_type, protein, position, reference, alternate, ambiguous_alternate, annotation, hgvs, variant_level) VALUES ('28262:G>GAACA', 'INSERTION', NULL, '28262', 'G', 'GAACA', FALSE, 'intergenic', NULL, 'nuc');

INSERT INTO lineage_variant_on(pango_lineage_id, variant_id) VALUES('P.1', '28262:G>GAACA');
