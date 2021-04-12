from typing import Tuple
from faker import Faker
from covigator.database.model import SampleEna, Sample, DataSource, JobEna, JobStatus, Log, CovigatorModule, Variant, \
    VariantObservation, VariantCooccurrence
from Bio.Alphabet.IUPAC import IUPACData


def get_mocked_variant(faker: Faker, chromosome=None, gene_name=None) -> Variant:
    return Variant(
        chromosome=chromosome if chromosome else faker.bothify(text="chr##"),
        position=faker.random_int(min=1, max=30000),
        reference=faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0],
        # TODO: reference and alternate could be equal!
        alternate=faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0],
        gene_name=gene_name,
        annotation=faker.random_choices(["non_synonymous", "inframe_deletion", "inframe_insertion"], length=1)[0]
    )


def get_mocked_variant_observation(sample: Sample, variant: Variant):
    return VariantObservation(
                    sample=sample.id,
                    source=sample.source,
                    chromosome=variant.chromosome,
                    position=variant.position,
                    reference=variant.reference,
                    alternate=variant.alternate
                )


def get_mocked_ena_sample(faker: Faker, job_status=JobStatus.COOCCURRENCE) -> Tuple[SampleEna, Sample, JobEna]:
    """
    Returns a triple of SampleEna, Sample and Job with the same sample identifier
    """
    identifier = faker.unique.uuid4()
    sample_ena = SampleEna(
        run_accession=identifier,
        first_created=faker.date_time(),
        country=faker.country(),
        fastq_ftp=faker.uri(),
        fastq_md5=faker.md5(),
        num_fastqs=1
    )
    sample = Sample(
        id=identifier,
        source=DataSource.ENA,
        ena_id=identifier
    )
    job = JobEna(
        run_accession=identifier,
        status=job_status
    )
    return sample_ena, sample, job


def get_mocked_log(faker: Faker, source: DataSource = None, module: CovigatorModule = None) -> Log:
    return Log(
        start=faker.date_time(),
        end=faker.date_time(),
        source=source if source else faker.random_choices((DataSource.ENA, DataSource.GISAID), length=1)[0],
        module=module if module else faker.random_choices((CovigatorModule.ACCESSOR, CovigatorModule.PROCESSOR),
                                                          length=1)[0],
        has_error=faker.boolean(),
        processed=faker.random_digit(),
        data={"included": faker.random_digit(),
              "excluded": {"this": faker.random_digit(), "that": faker.random_digit()}}
    )


def get_mocked_variant_cooccurrence(faker: Faker, variant_one: Variant, variant_two: Variant) -> VariantCooccurrence:
    return VariantCooccurrence(
        chromosome_one=variant_one.chromosome,
        position_one=variant_one.position,
        reference_one=variant_one.reference,
        alternate_one=variant_one.alternate,
        chromosome_two=variant_two.chromosome,
        position_two=variant_two.position,
        reference_two=variant_two.reference,
        alternate_two=variant_two.alternate,
        count=faker.random_int(min=1, max=10)
    )
