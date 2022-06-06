from itertools import combinations
from typing import Union, List

from sqlalchemy import and_
from sqlalchemy.orm import Session

from covigator import SYNONYMOUS_VARIANT
from covigator.database.model import JobStatus, VariantCooccurrence, GisaidVariantCooccurrence
from covigator.database.queries import Queries
from logzero import logger

BATCH_SIZE = 1000


class CooccurrenceMatrixLoader:

    def __init__(self, session: Session, source: str):
        self.session = session
        self.queries = Queries(session=self.session)
        self.source = source
        self.cache = {}
        self.variant_klazz = self.queries.get_variant_cooccurrence_klass(source=source)
        self.variant_cooccurrence_klazz = self.queries.get_variant_cooccurrence_klass(source=self.source)

    def load(self, maximum_length: int):

        # deletes the database before loading
        self.session.query(self.variant_cooccurrence_klazz).delete()

        # iterates over every sample in FINISHED status and computes the cooccurrence matrix
        sample_klass = self.queries.get_sample_klass(self.source)
        count_samples = self.queries.count_samples(source=self.source, cache=False)
        computed = 0
        query = self.session.query(sample_klass).filter(sample_klass.status == JobStatus.FINISHED)
        for sample in self.queries.windowed_query(query=query, column=sample_klass.run_accession, windowsize=BATCH_SIZE):
            self._compute_sample(sample.run_accession, maximum_length=maximum_length)
            computed += 1
            if computed % BATCH_SIZE == 0:
                # commits batches of 1000 samples
                self._commit_cache()
                logger.info('Processed cooccurrence over {}/{} ({} %) samples'.format(
                    computed, count_samples, round(float(computed) / count_samples * 100, 3)))

        # commits the last batch
        self._commit_cache()

        # once finished deletes the unique observations
        self.session.query(self.queries.get_variant_cooccurrence_klass(self.source)) \
            .filter(self.variant_cooccurrence_klazz.count == 1) \
            .delete()
        self.session.commit()

    def _get_variant_ids_by_sample(self, sample_id, source: str, maximum_length: int) -> List[str]:
        """
        Returns the variant ids of all mutations in a given sample after filtering out:
        mutations not overlapping any gene, synonymous mutations, long indels according to maximum_length parameter
        """
        klass = self.queries.get_variant_observation_klass(source=source)
        return self.session.query(klass.variant_id) \
            .filter(and_(klass.sample == sample_id,
                         klass.gene_name != None,
                         klass.annotation_highest_impact != SYNONYMOUS_VARIANT,
                         klass.length < maximum_length,
                         klass.length > -maximum_length)) \
            .order_by(klass.position, klass.reference, klass.alternate) \
            .all()

    def _unique_id(self, variant_id_one: str, variant_id_two: str):
        return "{}-{}".format(variant_id_one, variant_id_two)

    def _get_from_cache(self, variant_id_one: str, variant_id_two: str):
        return self.cache.get(self._unique_id(variant_id_one, variant_id_two), None)

    def _store_in_cache(self, variant_id_one: str, variant_id_two: str,
                        entry: Union[VariantCooccurrence, GisaidVariantCooccurrence]):
        self.cache[self._unique_id(variant_id_one, variant_id_two)] = entry

    def _commit_cache(self):
        self.session.add_all(list(self.cache.values()))
        self.session.commit()
        self.cache = {}

    def _increment_variant_cooccurrence(self, variant_id_one: str, variant_id_two: str):

        # NOTE: this method does not commit to DB due to performance reasons

        # first looks in the cache
        variant_cooccurrence = self._get_from_cache(variant_id_one=variant_id_one, variant_id_two=variant_id_two)
        if variant_cooccurrence is None:
            # if not in the cache looks in the DB
            variant_cooccurrence = self.session.query(self.variant_klazz) \
                .filter(and_(self.variant_klazz.variant_id_one == variant_id_one,
                             self.variant_klazz.variant_id_two == variant_id_two)) \
                .first()
        if variant_cooccurrence is None:
            # if not in the cache and not in the DB creates a new one
            variant_cooccurrence = self.variant_klazz(
                variant_id_one=variant_id_one,
                variant_id_two=variant_id_two,
                count=1)
        else:
            variant_cooccurrence.count = variant_cooccurrence.count + 1

        # stores the changes in the cache
        self._store_in_cache(variant_id_one=variant_id_one, variant_id_two=variant_id_two, entry=variant_cooccurrence)

    def _compute_sample(self, run_accession: str, maximum_length: int = 10):

        assert run_accession is not None or run_accession == "", "Missing sample identifier"

        sample_id = run_accession

        # the order by position is important to ensure we store only half the matrix and the same half of the matrix
        variant_ids = self._get_variant_ids_by_sample(sample_id, source=self.source, maximum_length=maximum_length)

        # process all pairwise combinations without repetitions including the diagoonal
        for (variant_id_one, variant_id_two) in list(combinations(variant_ids, 2)) + list(
                zip(variant_ids, variant_ids)):
            self._increment_variant_cooccurrence(variant_id_one, variant_id_two)
