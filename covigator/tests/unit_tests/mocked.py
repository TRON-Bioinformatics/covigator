import numpy as np
from itertools import combinations
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
    samples = mock_samples(faker=faker, session=session, num_samples=num_samples)
    for sample_ena, sample, job in samples:
        variants = [get_mocked_variant(faker=faker) for _ in range(10)]
        new_variants = list(filter(lambda x: x.variant_id not in existing_variants, variants))
        session.add_all(new_variants)
        session.commit()
        existing_variants.update([v.variant_id for v in variants])

        variants_observations = [get_mocked_variant_observation(faker=faker, variant=v, sample=sample)
                                 for v in variants]
        session.add_all(variants_observations)
        session.commit()


def mock_samples(faker, session: Session, num_samples=10, job_status=JobStatus.FINISHED):
    samples_ena = []
    samples = []
    jobs = []
    for _ in range(num_samples):
        sample_ena, sample, job = get_mocked_ena_sample(faker=faker, job_status=job_status)
        samples_ena.append(sample_ena)
        samples.append(sample)
        jobs.append(job)
    session.add_all(samples_ena)
    session.commit()
    session.add_all(samples)
    session.add_all(jobs)
    session.commit()
    return [(se, s, j) for se, s, j in zip(samples_ena, samples, jobs)]


def mock_cooccurrence_matrix(faker, session: Session):
    # add some variants belonging to two genes
    chromosome = "fixed_chromosome"
    gene_name = "S"
    variants = [get_mocked_variant(faker=faker, chromosome=chromosome, gene_name=gene_name) for _ in range(5)]
    other_gene_name = "N"
    other_variants = [get_mocked_variant(faker=faker, chromosome=chromosome, gene_name=other_gene_name) for _
                      in range(5)]
    session.add_all(variants + other_variants)
    session.commit()
    # adds some cooccurrences
    cooccurrences = []
    other_cooccurrences = []
    variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p): (v1, v2) for v1, v2 in
                          list(combinations(variants, 2))}
    other_variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p): (v1, v2) for v1, v2 in
                                list(combinations(other_variants, 2))}
    combined_variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p): (v1, v2) for v1, v2 in
                                   list(zip(variants, other_variants))}
    for variant in variants + other_variants:
        cooccurrences.append(get_mocked_variant_cooccurrence(faker, variant, variant))
    for (variant_one, variant_two) in [variants_to_sample.get(k) for k in
                                       np.random.choice(list(variants_to_sample.keys()), 5, replace=False)]:
        cooccurrences.append(get_mocked_variant_cooccurrence(faker, variant_one, variant_two))
    for (variant_one, variant_two) in [other_variants_to_sample.get(k) for k in
                                       np.random.choice(list(other_variants_to_sample.keys()), 5, replace=False)]:
        other_cooccurrences.append(get_mocked_variant_cooccurrence(faker, variant_one, variant_two))
    for (variant_one, variant_two) in [combined_variants_to_sample.get(k) for k in
                                       np.random.choice(list(combined_variants_to_sample.keys()), 5,
                                                        replace=False)]:
        other_cooccurrences.append(get_mocked_variant_cooccurrence(faker, variant_one, variant_two))
    session.add_all(cooccurrences + other_cooccurrences)
    session.commit()

    # add some samples to compute the frequency right
    mock_samples(faker=faker, session=session, num_samples=10)

    return other_variants, variants
