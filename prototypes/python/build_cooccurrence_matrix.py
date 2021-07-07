#!/bin/bash
from itertools import combinations
from covigator.database.database import Database
from covigator.database.model import SampleEna, VariantCooccurrence
from logzero import logger

from covigator.database.queries import Queries


def main():
    database = Database()
    session = database.get_database_session()
    queries = Queries(session=session)
    logger.info("Computing cooccurrence matrix")
    samples = session.query(SampleEna).all()
    processed = 0
    for sample in samples:
        sample_id = sample.run_accession
        processed += 1
        logger.info("Processing cooccurrent variants for sample {}/{}".format(processed, len(samples)))
        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variants = queries.get_variants_by_sample(sample_id)
        for (variant_one, variant_two) in combinations(variants, 2):
            variant_cooccurrence = queries.get_variant_cooccurrence(variant_one, variant_two)
            if variant_cooccurrence is None:
                variant_cooccurrence = VariantCooccurrence(
                    chromosome_one=variant_one.chromosome,
                    position_one=variant_one.position,
                    reference_one=variant_one.reference,
                    alternate_one=variant_one.alternate,
                    chromosome_two=variant_two.chromosome,
                    position_two=variant_two.position,
                    reference_two=variant_two.reference,
                    alternate_two=variant_two.alternate,
                    count=0
                )
                session.add(variant_cooccurrence)
            variant_cooccurrence.count += 1
        session.commit()


if __name__ == '__main__':
    main()
