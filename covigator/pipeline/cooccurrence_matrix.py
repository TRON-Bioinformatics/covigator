from itertools import combinations
from sqlalchemy.orm import Session
from covigator.database.queries import Queries


class CooccurrenceMatrixException(Exception):
    pass


class CooccurrenceMatrix:

    def compute(self, run_accession: str, source: str, session: Session):

        assert run_accession is not None or run_accession == "", "Missing sample identifier"
        assert session is not None, "Missing DB session"

        queries = Queries(session=session)
        sample_id = run_accession

        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variants = queries.get_variants_by_sample(sample_id, source=source)

        variant_cooccurrence_klass = queries.get_variant_cooccurrence_klass(source)

        # process all pairwise combinations without repetitions including the diagoonal
        for (variant_one, variant_two) in list(combinations(variants, 2)) + list(zip(variants, variants)):
            variant_cooccurrence = queries.get_variant_cooccurrence(variant_one, variant_two, source)
            if variant_cooccurrence is None:
                variant_cooccurrence = variant_cooccurrence_klass(
                    variant_id_one=variant_one.variant_id,
                    variant_id_two=variant_two.variant_id,
                    count=1
                )
                session.add(variant_cooccurrence)
            else:
                # NOTE: it is important to increase the counter like this to avoid race conditions
                # the increase happens in the database server and not in python
                # see https://stackoverflow.com/questions/2334824/how-to-increase-a-counter-in-sqlalchemy
                variant_cooccurrence.count = variant_cooccurrence_klass.count + 1
            session.commit()
