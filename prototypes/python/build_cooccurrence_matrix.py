#!/bin/bash
from itertools import combinations
from sqlalchemy import and_
from covigator.database import Database
from covigator.model import EnaRun, VariantObservation, VariantCooccurrence
from logzero import logger


def main():
    database = Database()
    session = database.get_database_session()
    logger.info("Computing cooccurrence matrix")
    for sample in session.query(EnaRun).all():
        sample_id = sample.run_accession
        logger.info("Processing cooccurrent variants for sample {}".format(sample_id))
        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variants = session.query(VariantObservation)\
            .filter(VariantObservation.sample == sample_id)\
            .order_by(VariantObservation.position)\
            .all()
        for (variant_one, variant_two) in combinations(variants, 2):
            variant_cooccurrence = session.query(VariantCooccurrence)\
                .filter(and_(VariantCooccurrence.chromosome_one == variant_one.chromosome,
                             VariantCooccurrence.position_one == variant_one.position,
                             VariantCooccurrence.reference_one == variant_one.reference,
                             VariantCooccurrence.alternate_one == variant_one.alternate,
                             VariantCooccurrence.chromosome_two == variant_two.chromosome,
                             VariantCooccurrence.position_two == variant_two.position,
                             VariantCooccurrence.reference_two == variant_two.reference,
                             VariantCooccurrence.alternate_two == variant_two.alternate))\
                .first()
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
