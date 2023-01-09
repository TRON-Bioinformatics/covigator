from typing import Union
from cyvcf2 import VCF, Variant
import os
import re
from sqlalchemy.exc import IntegrityError, InvalidRequestError
from sqlalchemy.orm import Session
from covigator import MISSENSE_VARIANT
from covigator.database.model import Variant as CovigatorVariant, VariantObservation, \
    SubclonalVariantObservation, SampleEna, DataSource, VariantType, \
    LowFrequencyVariantObservation, SubclonalVariant, LowFrequencyVariant, VariantCovid19Portal, \
    VariantObservationCovid19Portal, SampleCovid19Portal, LowQualityClonalVariant, LowQualityClonalVariantObservation
from covigator.database.queries import Queries
from covigator.exceptions import CovigatorNotSupportedVariant

MAX_INTRAHOST_SUBCLONAL = 5147

MAX_INTRAHOST_HIGH_FREQUENCY = 3

MIN_INTRAHOST_READ_COUNT = 50000

MIN_INTRAHOST_COVERED_BASES = 29000

MIN_INTRAHOST_LENGTH = 10

MIN_INTRAHOST_AC = 10

MIN_INTRAHOST_DP = 100


def is_eligible_intrahost_sample(sample, count_subclonal):
    return sample.covered_bases is not None and sample.covered_bases >= MIN_INTRAHOST_COVERED_BASES and \
        sample.read_count is not None and sample.read_count >= MIN_INTRAHOST_READ_COUNT and \
        count_subclonal <= MAX_INTRAHOST_SUBCLONAL


def is_eligible_intrahost_mutation(ac, dp, length, vaf):
    return vaf >= 0.02 and dp >= MIN_INTRAHOST_DP and ac >= MIN_INTRAHOST_AC and \
        length <= MIN_INTRAHOST_LENGTH


def _parse_additional_annotations(
        covigator_variant: Union[CovigatorVariant, SubclonalVariant, LowFrequencyVariant],
        vcf_variant: Variant):
    covigator_variant.cons_hmm_sars_cov_2 = vcf_variant.INFO.get("CONS_HMM_SARS_COV_2")
    covigator_variant.cons_hmm_sarbecovirus = vcf_variant.INFO.get("CONS_HMM_SARBECOVIRUS")
    covigator_variant.cons_hmm_vertebrate_cov = vcf_variant.INFO.get("CONS_HMM_VERTEBRATE_COV")
    covigator_variant.pfam_name = vcf_variant.INFO.get("PFAM_NAME")
    covigator_variant.pfam_description = vcf_variant.INFO.get("PFAM_DESCRIPTION")


def _parse_snpeff_annotations(
        covigator_variant: Union[CovigatorVariant, SubclonalVariant, LowFrequencyVariant],
        vcf_variant: Variant):
    ann = vcf_variant.INFO.get("ANN")
    if ann is not None:
        annotations = ann.split(",")
        annotation = annotations[0]  # NOTE: chooses arbitrarily the first annotation
        if "|" in annotation:
            values = annotation.split("|")
            covigator_variant.overlaps_multiple_genes = len(annotations) > 1
            covigator_variant.annotation = values[1].strip()
            covigator_variant.annotation_highest_impact = re.sub("&.*", "", covigator_variant.annotation)
            covigator_variant.annotation_impact = values[2].strip()
            covigator_variant.gene_name = values[3].strip()
            covigator_variant.gene_id = values[4].strip()
            covigator_variant.biotype = values[7].strip()
            covigator_variant.hgvs_c = values[9].strip()
            covigator_variant.hgvs_p = values[10].strip()
            covigator_variant.cdna_pos_length = values[11].strip()
            covigator_variant.cds_pos_length = values[12].strip()
            covigator_variant.aa_pos_length = values[13].strip()
            if covigator_variant.annotation == MISSENSE_VARIANT:
                hgvs_pattern = re.compile(r"^p\.([a-zA-Z]{1,3})([0-9]+)([a-zA-Z]{1,3})$")
                match = hgvs_pattern.match(covigator_variant.hgvs_p)
                if match:
                    covigator_variant.reference_amino_acid = match.group(1)
                    covigator_variant.alternate_amino_acid = match.group(3)
                    covigator_variant.position_amino_acid = int(match.group(2))


def _get_variant_type(reference, alternate):
    if len(reference) == 1 and len(alternate) == 1:
        return VariantType.SNV
    elif len(reference) == len(alternate):
        return VariantType.MNV
    elif len(reference) == 1 and len(alternate) > 1:
        return VariantType.INSERTION
    elif len(reference) > 1 and len(alternate) == 1:
        return VariantType.DELETION
    elif len(reference) > 1 and len(alternate) > 1:
        return VariantType.COMPLEX
    else:
        raise CovigatorNotSupportedVariant


def _get_variant_length(variant):
    return len(variant.ALT[0]) - len(variant.REF)


def _parse_variant_observation(
        variant: Variant, sample: Union[SampleEna, SampleCovid19Portal], covigator_variant: CovigatorVariant,
        klass):

    dp4 = variant.INFO.get("DP4")
    return klass(
        sample=sample.run_accession,
        variant_id=covigator_variant.variant_id,
        chromosome=variant.CHROM,
        position=variant.POS,
        reference=variant.REF,
        # NOTE: because variants are normalized we only expect one alternate
        alternate=variant.ALT[0],
        quality=variant.QUAL,
        filter=variant.FILTER,
        dp=variant.INFO.get("vafator_dp"),
        dp4_ref_forward=dp4[0] if dp4 else None,
        dp4_ref_reverse=dp4[1] if dp4 else None,
        dp4_alt_forward=dp4[2] if dp4 else None,
        dp4_alt_reverse=dp4[3] if dp4 else None,
        vaf=variant.INFO.get("vafator_af"),
        ac=variant.INFO.get("vafator_ac"),
        strand_bias=variant.INFO.get("SB"),
        # denormalized fields
        annotation=covigator_variant.annotation,
        annotation_impact=covigator_variant.annotation_impact,
        biotype=covigator_variant.biotype,
        annotation_highest_impact=covigator_variant.annotation_highest_impact,
        gene_name=covigator_variant.gene_name,
        hgvs_p=covigator_variant.hgvs_p,
        hgvs_c=covigator_variant.hgvs_c,
        date=sample.collection_date,
        variant_type=covigator_variant.variant_type,
        length=_get_variant_length(variant),
        reference_amino_acid=covigator_variant.reference_amino_acid,
        alternate_amino_acid=covigator_variant.alternate_amino_acid,
        position_amino_acid=covigator_variant.position_amino_acid,
        cons_hmm_sars_cov_2=covigator_variant.cons_hmm_sars_cov_2,
        cons_hmm_sarbecovirus=covigator_variant.cons_hmm_sarbecovirus,
        cons_hmm_vertebrate_cov=covigator_variant.cons_hmm_vertebrate_cov,
        pfam_name=covigator_variant.pfam_name,
        pfam_description=covigator_variant.pfam_description
    )


def _parse_variant(variant: Variant, klass) -> CovigatorVariant:
    parsed_variant = klass(
        chromosome=variant.CHROM,
        position=variant.POS,
        reference=variant.REF,
        # NOTE: because we only support haploid organisms we expect only one alternate,
        # TODO: support subclonal variants at some point
        alternate=variant.ALT[0],
        variant_type=_get_variant_type(reference=variant.REF, alternate=variant.ALT[0]),
        length=_get_variant_length(variant)
    )
    parsed_variant.variant_id = parsed_variant.get_variant_id()
    _parse_snpeff_annotations(covigator_variant=parsed_variant, vcf_variant=variant)
    _parse_additional_annotations(covigator_variant=parsed_variant, vcf_variant=variant)
    return parsed_variant


class VcfLoader:

    def load(self, vcf_file: str, run_accession: str, source: DataSource, session: Session):

        assert vcf_file is not None or vcf_file == "", "Missing VCF file provided to VcfLoader"
        assert os.path.exists(vcf_file) and os.path.isfile(vcf_file), "Non existing VCF file provided to VcfLoader"
        assert run_accession is not None or run_accession == "", "Missing sample"
        assert session is not None, "Missing DB session"

        observed_variants = []
        subclonal_observed_variants = []
        low_frequency_observed_variants = []
        specific_sample = Queries(session=session).find_sample_by_accession(
            run_accession=run_accession, source=source)
        assert specific_sample is not None, "Cannot find sample in database"

        # reads whole VCF in memory to count variants
        variants = [v for v in VCF(vcf_file)]

        # counts variants
        specific_sample.count_snvs = len([v for v in variants if v.FILTER is None and len(v.REF) == 1 and len(v.ALT[0]) == 1])
        specific_sample.count_deletions = len([v for v in variants if v.FILTER is None and len(v.REF) == 1 and len(v.ALT[0]) > 1])
        specific_sample.count_insertions = len([v for v in variants if v.FILTER is None and len(v.REF) > 1 and len(v.ALT[0]) == 1])

        count_high_frequency_non_clonal = 0
        count_subclonal = 0
        if source == DataSource.ENA:
            specific_sample.count_subclonal_snvs = len([v for v in variants if v.FILTER == "SUBCLONAL"
                                        and len(v.REF) == 1 and len(v.ALT[0]) == 1])
            specific_sample.count_subclonal_deletions = len([v for v in variants if v.FILTER == "SUBCLONAL"
                                             and len(v.REF) > 1 and len(v.ALT[0]) == 1])
            specific_sample.count_subclonal_insertions = len([v for v in variants if v.FILTER == "SUBCLONAL"
                                              and len(v.REF) == 1 and len(v.ALT[0]) > 1])
            specific_sample.count_low_frequency_snvs = len([v for v in variants if v.FILTER  == "LOW_FREQUENCY"
                                        and len(v.REF) == 1 and len(v.ALT[0]) == 1])
            specific_sample.count_low_frequency_deletions = len([v for v in variants if v.FILTER == "LOW_FREQUENCY"
                                             and len(v.REF) > 1 and len(v.ALT[0]) == 1])
            specific_sample.count_low_frequency_insertions = len([v for v in variants if v.FILTER == "LOW_FREQUENCY"
                                              and len(v.REF) == 1 and len(v.ALT[0]) > 1])

            count_subclonal = len([v for v in variants if v.INFO.get("vafator_af", 0.0) < 0.8])
            specific_sample.intrahost_filter = not is_eligible_intrahost_sample(
                sample=specific_sample, count_subclonal=count_subclonal)

            count_high_frequency_non_clonal = len([v for v in variants if 0.4 <= v.INFO.get("vafator_af", 0.0) < 0.8])
            specific_sample.potential_coinfection = count_high_frequency_non_clonal > MAX_INTRAHOST_HIGH_FREQUENCY


        variant: Variant
        for variant in variants:
            # NOTE: it accepts only variant flagged as PASS, LOW_FREQUENCY, SUBCLONAL
            # but this classification is irrelevant as it then classifies mutations based on VAF
            # furthermore, the classification from the pipeline may not consistent with the VAF
            if variant.FILTER is None or variant.FILTER in ["LOW_FREQUENCY", "SUBCLONAL", "LOW_QUALITY_CLONAL"]:
                if source == DataSource.COVID19_PORTAL:
                    v = _parse_variant(variant, VariantCovid19Portal)
                    observed_variants.append(
                        _parse_variant_observation(
                            variant, specific_sample, v, VariantObservationCovid19Portal))
                    session.add(v)
                else:
                    vaf = variant.INFO.get("vafator_af", 0.0)
                    dp = variant.INFO.get("vafator_dp", 0)
                    ac = variant.INFO.get("vafator_ac", 0)
                    length = abs(_get_variant_length(variant))
                    if vaf >= 0.8:
                        # only stores clonal high quality variants in this table
                        v = _parse_variant(variant, CovigatorVariant)
                        observed_variants.append(
                            _parse_variant_observation(
                                variant, specific_sample, v, VariantObservation))
                        session.add(v)
                    elif vaf >= 0.5:  # and < 0.8
                        # stores clonal low quality variants in this table
                        v = _parse_variant(variant, LowQualityClonalVariant)
                        observed_variants.append(
                            _parse_variant_observation(
                                variant, specific_sample, v, LowQualityClonalVariantObservation))
                        session.add(v)
                    elif is_eligible_intrahost_mutation(ac=ac, dp=dp, length=length, vaf=vaf) and \
                            not specific_sample.intrahost_filter and \
                            not specific_sample.potential_coinfection:  # and vaf < 0.5
                        # stores high quality intrahost subclonal variants in this table
                        v = _parse_variant(variant, SubclonalVariant)
                        subclonal_observed_variants.append(
                            _parse_variant_observation(
                                variant, specific_sample, v, SubclonalVariantObservation))
                        session.add(v)
                    else:
                        # stores low frequency and/or low quality subclonal variants in this table (we do nothing with these variants)
                        v = _parse_variant(variant, LowFrequencyVariant)
                        low_frequency_observed_variants.append(
                            _parse_variant_observation(
                                variant, specific_sample, v, LowFrequencyVariantObservation))
                        session.add(v)

                try:
                    session.commit()
                except (IntegrityError, InvalidRequestError):
                    # do nothing the variant was just added by another process between merge and commit
                    session.rollback()
        session.add_all(observed_variants)
        session.add_all(subclonal_observed_variants)
        session.add_all(low_frequency_observed_variants)
        # NOTE: commit will happen afterwards when the job status is updated
