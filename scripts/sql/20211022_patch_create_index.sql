

CREATE INDEX variant_observation_test_index_annotation_position ON public.variant_observation_test USING btree (annotation_highest_impact, "position");
CREATE INDEX variant_observation_test_index_sample ON public.variant_observation_test USING btree (sample);