from itertools import combinations
from sqlalchemy.orm import Session
from covigator.database.model import Sample, VariantCooccurrence
from logzero import logger
from covigator.database.queries import Queries
from sqlalchemy.exc import IntegrityError


class CooccurrenceMatrixException(Exception):
    pass


class CooccurrenceMatrix:

    def compute(self, sample: Sample, session: Session):

        assert sample is not None, "Missing sample"
        assert sample.id is not None or sample.id == "", "Missing sample identifier"
        assert session is not None, "Missing DB session"

        queries = Queries(session=session)
        sample_id = sample.id
        logger.info("Processing cooccurrent variants for sample {}".format(sample_id))

        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variants = queries.get_variants_by_sample(sample_id)
        failed_variants = []

        # process all pairwise combinations without repetitions including the diagoonal
        for (variant_one, variant_two) in list(combinations(variants, 2)) + list(zip(variants, variants)):
            try:
                variant_cooccurrence = queries.get_variant_cooccurrence(variant_one, variant_two)
                if variant_cooccurrence is None:
                    variant_cooccurrence = VariantCooccurrence(
                        variant_id_one=variant_one.variant_id,
                        variant_id_two=variant_two.variant_id,
                        count=1
                    )
                    session.add(variant_cooccurrence)
                    session.commit()
                else:
                    # NOTE: it is important to increase the counter like this to avoid race conditions
                    # the increase happens in the database server and not in python
                    # see https://stackoverflow.com/questions/2334824/how-to-increase-a-counter-in-sqlalchemy
                    variant_cooccurrence.count = VariantCooccurrence.count + 1
            except IntegrityError:
                session.rollback()
                failed_variants.append((variant_one, variant_two))

        # tries again the failed variants as these are expected to be there now
        for (variant_one, variant_two) in failed_variants:
            variant_cooccurrence = queries.get_variant_cooccurrence(variant_one, variant_two)
            if variant_cooccurrence is None:
                raise CooccurrenceMatrixException("Some cooccurrent variants failed to be persisted twice")
            variant_cooccurrence.count = VariantCooccurrence.count + 1
