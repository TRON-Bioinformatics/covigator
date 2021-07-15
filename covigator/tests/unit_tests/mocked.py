from typing import Tuple
from faker import Faker
from sqlalchemy.orm import Session

from covigator.database.model import SampleEna, Sample, DataSource, JobEna, JobStatus, Log, CovigatorModule, Variant, \
    VariantObservation, VariantCooccurrence, VariantType
from Bio.Alphabet.IUPAC import IUPACData


def get_mocked_variant(faker: Faker, chromosome=None, gene_name=None) -> Variant:
    variant = Variant(
        chromosome=chromosome if chromosome else faker.bothify(text="chr##"),
        position=faker.random_int(min=1, max=30000),
        reference=faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0],
        # TODO: reference and alternate could be equal!
        alternate=faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0],
        variant_type=VariantType.SNV,
        gene_name=gene_name,
        hgvs_p="p.{}{}{}".format(
            faker.random_choices(list(IUPACData.protein_letters_1to3.values()), length=1)[0],
            faker.random_int(min=1, max=500),
            faker.random_choices(list(IUPACData.protein_letters_1to3.values()), length=1)[0]
        ),
        annotation=faker.random_choices(["non_synonymous", "inframe_deletion", "inframe_insertion"], length=1)[0]
    )
    variant.variant_id = variant.get_variant_id()
    return variant


def get_mocked_variant_observation(sample: Sample, variant: Variant, faker=Faker()):
    return VariantObservation(
        sample=sample.id if sample else faker.unique.uuid4(),
        source=sample.source if sample else faker.random_choices((DataSource.ENA, DataSource.GISAID)),
        variant_id=variant.variant_id,
        chromosome=variant.chromosome,
        position=variant.position,
        reference=variant.reference,
        alternate=variant.alternate,
        variant_type=VariantType.SNV)


def get_mocked_ena_sample(faker: Faker, job_status=JobStatus.FINISHED) -> Tuple[SampleEna, Sample, JobEna]:
    """
    Returns a triple of SampleEna, Sample and Job with the same sample identifier
    """
    identifier = faker.unique.uuid4()
    sample_ena = SampleEna(
        run_accession=identifier,
        collection_date=faker.date_time(),
        country=faker.country(),
        fastq_ftp=faker.uri(),
        fastq_md5=faker.md5(),
        num_fastqs=1,
        finished=job_status == JobStatus.FINISHED
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
    if variant_one.position <= variant_two.position:
        cooccurrence = VariantCooccurrence(
            variant_id_one=variant_one.variant_id,
            variant_id_two=variant_two.variant_id,
            count=faker.random_int(min=1, max=10)
        )
    else:
        cooccurrence = VariantCooccurrence(
            variant_id_two=variant_one.variant_id,
            variant_id_one=variant_two.variant_id,
            count=faker.random_int(min=1, max=10)
        )
    return cooccurrence


def mock_samples_and_variants(faker, session: Session, num_samples=10):
    existing_variants = set()
    for _ in range(num_samples):
        sample_ena, sample, job = get_mocked_ena_sample(faker=faker)
        session.add_all([sample_ena, sample, job])
        variants = [get_mocked_variant(faker=faker) for _ in range(10)]
        session.add_all(list(filter(lambda x: x in existing_variants, variants)))
        existing_variants.update([v.variant_id for v in variants])
        session.commit()
        variants_observations = [get_mocked_variant_observation(faker=faker, variant=v, sample=sample)
                                 for v in variants]
        session.add_all(variants_observations)
    session.commit()
