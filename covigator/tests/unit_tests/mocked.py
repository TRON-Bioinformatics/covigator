import numpy as np
from itertools import combinations
from typing import Tuple, Union
from faker import Faker
from sqlalchemy.orm import Session

from covigator import MISSENSE_VARIANT, SYNONYMOUS_VARIANT, INFRAME_INSERTION, INFRAME_DELETION
from covigator.database.model import SampleEna, DataSource, JobStatus, Log, CovigatorModule, Variant, \
    VariantCooccurrence, VariantType, SampleGisaid, GisaidVariant
from Bio.Alphabet.IUPAC import IUPACData

from covigator.database.queries import Queries

MOCKED_GENES = ["S", "N", "E", "M", "ORF3a", "ORF1ab", "ORF7b", "ORF10", "ORF6", "ORF8", "ORF7a"]
MOCKED_DOMAINS = ["CoV_S2", "bCoV_S1_RBD", "bCoV_S1_N", "CoV_S1_C"]
MOCKED_LINEAGES = ["A", "B", "C", "D", "E", "F", ""]


def get_mocked_variant(faker: Faker, chromosome=None, gene_name=None, source=DataSource.ENA.name, session: Session = None) -> Variant:

    if gene_name is None:
        gene_name = faker.random_choices(MOCKED_GENES, length=1)[0]
    domain_name = faker.random_choices(MOCKED_DOMAINS, length=1)[0]
    annotation = faker.random_choices(
            [MISSENSE_VARIANT, INFRAME_DELETION, INFRAME_INSERTION], length=1)[0] #SYNONYMOUS_VARIANT

    if session:
        queries = Queries(session=session)
        gene = queries.get_gene(gene_name=gene_name)
        start = gene.start
        end = gene.end
    else:
        start = 1
        end = 30000

    klass = Queries.get_variant_klass(source)
    reference = faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0]
    alternate = faker.random_choices(list(IUPACData.unambiguous_dna_letters), length=1)[0]
    variant = klass(
        chromosome=chromosome if chromosome else faker.bothify(text="chr##"),
        position=faker.random_int(min=start, max=end),
        reference=reference,
        # TODO: reference and alternate could be equal!
        alternate=alternate,
        variant_type=VariantType.SNV,
        gene_name=gene_name,
        hgvs_p="p.{}{}{}".format(
            faker.random_choices(list(IUPACData.protein_letters_1to3.values()), length=1)[0],
            faker.random_int(min=1, max=500),
            faker.random_choices(list(IUPACData.protein_letters_1to3.values()), length=1)[0]
        ),
        annotation=annotation,
        annotation_highest_impact=annotation,
        pfam_name=domain_name,
        length=len(reference) - len(alternate)
    )
    variant.variant_id = variant.get_variant_id()
    return variant


def get_mocked_variant_observation(
        sample: Union[SampleEna, SampleGisaid], variant: Union[Variant or GisaidVariant], faker=Faker()):

    klass = Queries.get_variant_observation_klass(
        DataSource.ENA.name if isinstance(sample, SampleEna) else DataSource.GISAID.name)
    return klass(
        sample=sample.run_accession if sample else faker.unique.uuid4(),
        variant_id=variant.variant_id,
        chromosome=variant.chromosome,
        position=variant.position,
        reference=variant.reference,
        alternate=variant.alternate,
        variant_type=VariantType.SNV,
        annotation=variant.annotation,
        annotation_highest_impact=variant.annotation_highest_impact,
        gene_name=variant.gene_name,
        pfam_name=variant.pfam_name,
        date=faker.date_time(),
        hgvs_p=variant.hgvs_p,
        length=len(variant.reference) - len(variant.alternate)
    )


def get_mocked_sample(faker: Faker, source: DataSource, job_status=JobStatus.FINISHED) -> Union[SampleEna, SampleGisaid]:
    identifier = faker.unique.uuid4()
    if source == DataSource.ENA:
        sample = SampleEna(
            run_accession=identifier,
            collection_date=faker.date_time(),
            country=faker.country(),
            fastq_ftp=faker.uri(),
            fastq_md5=faker.md5(),
            num_fastqs=1,
            status=job_status,
            pangolin_lineage=faker.random_choices(MOCKED_LINEAGES, length=1)[0]
        )
    elif source == DataSource.GISAID:
        sample = SampleGisaid(
            run_accession=identifier,
            collection_date=faker.date_time(),
            country=faker.country(),
            status=job_status,
            pangolin_lineage=faker.random_choices(MOCKED_LINEAGES, length=1)[0]
        )
    else:
        raise ValueError("Bad data source")
    return sample


def get_mocked_ena_sample(faker: Faker, job_status=JobStatus.FINISHED) -> SampleEna:
    return get_mocked_sample(faker=faker, source=DataSource.ENA, job_status=job_status)


def get_mocked_gisaid_sample(faker: Faker, job_status=JobStatus.FINISHED) -> SampleGisaid:
    return get_mocked_sample(faker=faker, source=DataSource.GISAID, job_status=job_status)


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
    existing_variants = {DataSource.ENA.name: set(), DataSource.GISAID.name: set()}
    samples = mock_samples(faker=faker, session=session, num_samples=num_samples)
    for sample in samples:
        source = DataSource.ENA if isinstance(sample, SampleEna) else DataSource.GISAID
        variants = [get_mocked_variant(faker=faker, source=source.name, session=session) for _ in range(10)]
        # this aims at removing potentially repeated variants
        variants_dict = {v.variant_id: v for v in variants}
        variants = variants_dict.values()
        new_variants = list(filter(lambda x: x.variant_id not in existing_variants.get(source.name), variants))
        session.add_all(new_variants)
        session.commit()
        existing_variants.get(source.name).update([v.variant_id for v in variants])

        variants_observations = [get_mocked_variant_observation(faker=faker, variant=v, sample=sample)
                                 for v in variants]
        session.add_all(variants_observations)
        session.commit()


def mock_samples(faker, session: Session, num_samples=10, job_status=JobStatus.FINISHED, source=None):
    samples = []
    for _ in range(num_samples):
        if source is not None:
            selected_source = source
        else:
            selected_source = faker.random_choices([DataSource.ENA.name, DataSource.GISAID.name], length=1)[0]
        if selected_source == DataSource.ENA.name:
            sample= get_mocked_ena_sample(faker=faker, job_status=job_status)
        else:   # GISAID
            sample = get_mocked_gisaid_sample(faker=faker, job_status=job_status)
        samples.append(sample)

    session.add_all(samples)
    session.commit()
    return samples


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
