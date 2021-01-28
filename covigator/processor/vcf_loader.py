from cyvcf2 import VCF, Variant
import os
from sqlalchemy.orm import Session
from covigator.model import Variant as CovigatorVariant, VariantObservation


class VcfLoader:

    def load(self, vcf_file: str, sample: str, session: Session):

        assert vcf_file is not None or vcf_file == "", "Missing VCF file provided to VcfLoader"
        assert os.path.exists(vcf_file) and os.path.isfile(vcf_file), "Non existing VCF file provided to VcfLoader"
        assert sample is not None or sample == "", "Missing sample"
        assert session is not None, "Missing DB session"

        observed_variants = []
        variant: Variant
        for variant in VCF(vcf_file):

            covigator_variant = self._parse_variant(variant)
            if covigator_variant:
                # NOTE: merge checks for existence adds or updates it if required
                # this variant is not part of the rollback if something else fails
                session.merge(covigator_variant)
            session.commit()
            observed_variants.append(self._parse_variant_observation(variant, sample))

        session.add_all(observed_variants)  # NOTE: commit will happen afterwards when the job status is updated

    def _parse_variant(self, variant: Variant) -> CovigatorVariant:
        parsed_variant = CovigatorVariant(
                    chromosome=variant.CHROM,
                    position=variant.POS,
                    reference=variant.REF,
                    # NOTE: because we only support haploid organisms we expect only one alternate,
                    # TODO: support subclonal variants at some point
                    alternate=variant.ALT[0])
        ann = variant.INFO.get("ANN")
        if ann is not None:
            annotations = ann.split(",")
            annotation = annotations[0]     # NOTE: chooses arbitrarily the first annotation
            if "|" in annotation:
                values = annotation.split("|")
                parsed_variant.overlaps_multiple_genes=len(annotations) > 1
                parsed_variant.annotation=values[1].strip()
                parsed_variant.annotation_impact=values[2].strip()
                parsed_variant.gene_name=values[3].strip()
                parsed_variant.gene_id=values[4].strip()
                parsed_variant.biotype=values[7].strip()
                parsed_variant.hgvs_c=values[9].strip()
                parsed_variant.hgvs_p=values[10].strip()
                parsed_variant.cdna_pos_length=values[11].strip()
                parsed_variant.cds_pos_length=values[12].strip()
                parsed_variant.aa_pos_length=values[13].strip()
        return parsed_variant

    def _parse_variant_observation(self, variant: Variant, sample: str) -> VariantObservation:

        dp4 = variant.INFO.get("DP4")
        return VariantObservation(
            sample=sample,
            chromosome=variant.CHROM,
            position=variant.POS,
            reference=variant.REF,
            # NOTE: because we only support haploid organisms we expect only one alternate,
            # TODO: support sublconal variants at some point
            alternate=variant.ALT[0],
            quality=variant.QUAL,
            filter=variant.FILTER,
            dp=variant.INFO.get("DP"),
            dp4_ref_forward=dp4[0] if dp4 else None,
            dp4_ref_reverse=dp4[1] if dp4 else None,
            dp4_alt_forward=dp4[2] if dp4 else None,
            dp4_alt_reverse=dp4[3] if dp4 else None,
            mapping_quality=variant.INFO.get("MQ"),
            mapping_quality_zero_fraction=variant.INFO.get("MQ0F")
        )
