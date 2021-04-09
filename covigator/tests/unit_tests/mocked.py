from typing import Tuple

from faker import Faker

from covigator.database.model import SampleEna, Sample, DataSource, JobEna, JobStatus, Log, CovigatorModule, Variant, \
    VariantObservation


def get_mocked_variant(faker: Faker) -> Variant:
    return Variant(
        chromosome=faker.bothify(text="chr##"),
        position=faker.random_int(min=1, max=30000),
        reference=faker.random_choices(["A", "C", "G", "T"], length=1)[0],
        # TODO: reference and alternate could be equal!
        alternate=faker.random_choices(["A", "C", "G", "T"], length=1)[0]
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


def get_mocked_ena_sample(faker: Faker, job_status=JobStatus.LOADED) -> Tuple[SampleEna, Sample, JobEna]:
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
