from typing import Union
from cyvcf2 import VCF, Variant
import os
import re
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session
from covigator.database.model import Variant as CovigatorVariant, VariantObservation, Sample, \
    SubclonalVariantObservation, SampleEna, SampleGisaid, DataSource, VariantType
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorNotSupportedVariant, CovigatorExcludedSampleTooManyMutations


class VcfLoader:

    def load(self, vcf_file: str, sample: Sample, session: Session,
             max_snvs=None, max_deletions=None, max_insertions=None):

        assert vcf_file is not None or vcf_file == "", "Missing VCF file provided to VcfLoader"
        assert os.path.exists(vcf_file) and os.path.isfile(vcf_file), "Non existing VCF file provided to VcfLoader"
        assert sample.id is not None or sample.id == "", "Missing sample"
        assert session is not None, "Missing DB session"

        observed_variants = []
        subclonal_observed_variants = []
        specific_sample = Queries(session=session).find_sample_by_accession(
            run_accession=sample.id, source=sample.source)
        assert specific_sample is not None, "Cannot find sample in database"

        # reads whole VCF in memory to count variants
        variants = [v for v in VCF(vcf_file)]

        # counts variants
        count_snvs = len([v for v in variants if v.FILTER is None and len(v.REF) == 1 and len(v.ALT[0]) == 1])
        count_insertions = len([v for v in variants if v.FILTER is None and len(v.REF) == 1 and len(v.ALT[0]) > 1])
        count_deletions = len([v for v in variants if v.FILTER is None and len(v.REF) > 1 and len(v.ALT[0]) == 1])
        too_many_snvs = max_snvs is not None and count_snvs > max_snvs
        too_many_insertions = max_insertions is not None and count_insertions > max_insertions
        too_many_deletions = max_deletions is not None and count_deletions > max_deletions
        if too_many_snvs or too_many_deletions or too_many_insertions:
            raise CovigatorExcludedSampleTooManyMutations(
                "Too many variants: SNVs={}, insertions={}, deletions={}".format(
                    count_snvs, count_insertions, count_deletions))

        specific_sample.count_snvs = count_snvs
        specific_sample.count_deletions = count_deletions
        specific_sample.count_insertions = count_insertions

        if sample.source == DataSource.ENA:
            count_subclonal_snvs = len([v for v in variants if v.FILTER in ["LOW_FREQUENCY", "SUBCLONAL"]
                                        and len(v.REF) == 1 and len(v.ALT[0]) == 1])
            count_subclonal_deletions = len([v for v in variants if v.FILTER in ["LOW_FREQUENCY", "SUBCLONAL"]
                                             and len(v.REF) > 1 and len(v.ALT[0]) == 1])
            count_subclonal_insertions = len([v for v in variants if v.FILTER in ["LOW_FREQUENCY", "SUBCLONAL"]
                                              and len(v.REF) == 1 and len(v.ALT[0]) > 1])
            specific_sample.count_subclonal_snvs = count_subclonal_snvs
            specific_sample.count_subclonal_deletions = count_subclonal_deletions
            specific_sample.count_subclonal_insertions = count_subclonal_insertions

        variant: Variant
        for variant in variants:
            if variant.FILTER is None or variant.FILTER in ["LOW_FREQUENCY", "SUBCLONAL"]:
                covigator_variant = self._parse_variant(variant)
                if covigator_variant:
                    # NOTE: merge checks for existence adds or updates it if required
                    # this variant is not part of the rollback if something else fails
                    session.merge(covigator_variant)
                    try:
                        session.commit()
                    except IntegrityError:
                        # do nothing the variant was just added by another process between merge and commit
                        pass
                if variant.FILTER is None:
                    # only stores clonal high quality variants in this table
                    observed_variants.append(
                        self._parse_variant_observation(variant, specific_sample, sample.source, covigator_variant, VariantObservation))
                elif variant.FILTER in ["LOW_FREQUENCY", "SUBCLONAL"]:
                    subclonal_observed_variants.append(
                        self._parse_variant_observation(variant, specific_sample, sample.source, covigator_variant, SubclonalVariantObservation))
        session.add_all(observed_variants)
        session.add_all(subclonal_observed_variants)
        # NOTE: commit will happen afterwards when the job status is updated

    def _parse_variant(self, variant: Variant) -> CovigatorVariant:
        parsed_variant = CovigatorVariant(
            chromosome=variant.CHROM,
            position=variant.POS,
            reference=variant.REF,
            # NOTE: because we only support haploid organisms we expect only one alternate,
            # TODO: support subclonal variants at some point
            alternate=variant.ALT[0],
            variant_type=self._get_variant_type(reference=variant.REF, alternate=variant.ALT[0]),
            length=self._get_variant_length(variant)
        )
        parsed_variant.variant_id = parsed_variant.get_variant_id()
        ann = variant.INFO.get("ANN")
        if ann is not None:
            annotations = ann.split(",")
            annotation = annotations[0]     # NOTE: chooses arbitrarily the first annotation
            if "|" in annotation:
                values = annotation.split("|")
                parsed_variant.overlaps_multiple_genes = len(annotations) > 1
                parsed_variant.annotation = values[1].strip()
                parsed_variant.annotation_highest_impact = re.sub("&.*", "", parsed_variant.annotation)
                parsed_variant.annotation_impact = values[2].strip()
                parsed_variant.gene_name = values[3].strip()
                parsed_variant.gene_id = values[4].strip()
                parsed_variant.biotype = values[7].strip()
                parsed_variant.hgvs_c = values[9].strip()
                parsed_variant.hgvs_p = values[10].strip()
                parsed_variant.cdna_pos_length = values[11].strip()
                parsed_variant.cds_pos_length = values[12].strip()
                parsed_variant.aa_pos_length = values[13].strip()
                if parsed_variant.annotation == "missense_variant":
                    hgvs_pattern = re.compile(r"^p\.([a-zA-Z]{1,3})([0-9]+)([a-zA-Z]{1,3})$")
                    match = hgvs_pattern.match(parsed_variant.hgvs_p)
                    if match:
                        parsed_variant.reference_amino_acid = match.group(1)
                        parsed_variant.alternate_amino_acid = match.group(3)
                        parsed_variant.position_amino_acid = int(match.group(2))
        return parsed_variant

    def _parse_variant_observation(
            self, variant: Variant, sample: Union[SampleEna, SampleGisaid], source: DataSource, covigator_variant: CovigatorVariant, klass):

        dp4 = variant.INFO.get("DP4")
        return klass(
            sample=sample.run_accession,
            source=source,
            variant_id=covigator_variant.variant_id,
            chromosome=variant.CHROM,
            position=variant.POS,
            reference=variant.REF,
            # NOTE: because variants are normalized we only expect one alternate
            alternate=variant.ALT[0],
            quality=variant.QUAL,
            filter=variant.FILTER,
            dp=variant.INFO.get("DP"),
            dp4_ref_forward=dp4[0] if dp4 else None,
            dp4_ref_reverse=dp4[1] if dp4 else None,
            dp4_alt_forward=dp4[2] if dp4 else None,
            dp4_alt_reverse=dp4[3] if dp4 else None,
            vaf=variant.INFO.get("AF"),
            strand_bias=variant.INFO.get("SB"),
            # denormalized fields
            annotation=covigator_variant.annotation,
            gene_name=covigator_variant.gene_name,
            hgvs_p=covigator_variant.hgvs_p,
            hgvs_c=covigator_variant.hgvs_c,
            date=sample.first_created if source == DataSource.ENA else sample.date,
            variant_type=self._get_variant_type(reference=variant.REF, alternate=variant.ALT[0]),
            length=self._get_variant_length(variant),
            reference_amino_acid=covigator_variant.reference_amino_acid,
            alternate_amino_acid=covigator_variant.alternate_amino_acid,
            position_amino_acid=covigator_variant.position_amino_acid
        )

    def _get_variant_type(self, reference, alternate):
        if len(reference) and len(alternate):
            return VariantType.SNV
        elif len(reference) == 1 and len(alternate) > 1:
            return VariantType.INSERTION
        elif len(reference) > 1 and len(alternate) == 1:
            return VariantType.DELETION
        else:
            raise CovigatorNotSupportedVariant

    def _get_variant_length(self, variant):
        return len(variant.ALT[0]) - len(variant.REF)
