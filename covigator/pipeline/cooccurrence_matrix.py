from itertools import combinations
from sqlalchemy.orm import Session
from covigator.database.queries import Queries


class CooccurrenceMatrix:

    def compute(self, run_accession: str, source: str, session: Session, maximum_length: int = 10):

        assert run_accession is not None or run_accession == "", "Missing sample identifier"
        assert session is not None, "Missing DB session"

        queries = Queries(session=session)
        sample_id = run_accession

        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variant_ids = queries.get_variant_ids_by_sample(sample_id, source=source, maximum_length=maximum_length)

        # process all pairwise combinations without repetitions including the diagoonal
        new_coocurrences = []
        for (variant_id_one, variant_id_two) in list(combinations(variant_ids, 2)) + list(zip(variant_ids, variant_ids)):
            vc = queries.increment_variant_cooccurrence(variant_id_one, variant_id_two, source)
            if vc is not None:
                new_coocurrences.append(vc)

        session.add_all(new_coocurrences)
        session.commit()
