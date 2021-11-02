

CREATE INDEX variant_observation_v19_index_annotation_position ON public.variant_observation_v19 USING btree (annotation_highest_impact, "position");
CREATE INDEX variant_observation_v19_index_sample ON public.variant_observation_v19 USING btree (sample);
CREATE INDEX variant_observation_v19_index_position ON public.variant_observation_v19 USING btree (position);
CREATE INDEX variant_observation_v19_index_annotation_source ON public.variant_observation_v19 USING btree (annotation_highest_impact, "source");
CREATE INDEX variant_observation_v19_index_variant_id_source ON public.variant_observation_v19 USING btree (variant_id, "source");


CREATE INDEX subclonal_variant_observation_v19_index_position ON public.subclonal_variant_observation_v19 USING btree (position);
CREATE INDEX subclonal_variant_observation_v19_index_annotation_vaf ON public.subclonal_variant_observation_v19 USING btree (annotation_highest_impact, vaf);
CREATE INDEX subclonal_variant_observation_v19_index_vaf ON public.subclonal_variant_observation_v19 USING btree (vaf);