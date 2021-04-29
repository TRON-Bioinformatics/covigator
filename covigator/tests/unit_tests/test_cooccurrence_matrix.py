from unittest import TestCase

from faker import Faker

from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.database.model import Sample, VariantCooccurrence
from covigator.pipeline.cooccurrence_matrix import CooccurrenceMatrix
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample, get_mocked_variant, get_mocked_variant_observation


class CooccurrenceMatrixTests(TestCase):

    NUM_UNIQUE_VARIANTS = 10
    NUM_SAMPLES = 10
    NUM_VARIANT_OBSERVATIONS_PER_SAMPLE = 5

    def setUp(self) -> None:
        self.config = Configuration()
        self.session = Database(test=True, config=self.config).get_database_session()
        faker = Faker()
        # mocks some unique variants
        self.mocked_variants = [get_mocked_variant(faker) for _ in range(self.NUM_UNIQUE_VARIANTS)]
        self.session.add_all(self.mocked_variants)
        self.session.commit()

        variant_observations = []
        self.samples = []
        for _ in range(self.NUM_SAMPLES):
            sample_ena, sample, job = get_mocked_ena_sample(faker)
            self.session.add(sample)
            self.session.add(sample_ena)
            self.session.add(job)
            self.session.commit()
            self.samples.append(sample)
            # mocks observed variants
            for vo in faker.random_elements(
                    self.mocked_variants, length=self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE, unique=True):
                variant_observations.append(get_mocked_variant_observation(sample, vo))
            self.session.add_all(variant_observations)
            self.session.commit()

    def test_missing_sample(self):
        self.assertRaises(
            AssertionError,
            CooccurrenceMatrix().compute,
            None,
            self.session)

    def test_missing_session(self):
        self.assertRaises(
            AssertionError,
            CooccurrenceMatrix().compute,
            Sample(id="12345"),
            None)

    def test_non_existing_sample_does_not_add_new_entries(self):
        count = self.session.query(VariantCooccurrence).count()
        self.assertEqual(count, 0)
        CooccurrenceMatrix().compute(Sample(id="12345"), self.session)
        self.session.commit()
        count = self.session.query(VariantCooccurrence).count()
        self.assertEqual(count, 0)

    def test_one_existing_sample(self):
        count = self.session.query(VariantCooccurrence).count()
        self.assertEqual(count, 0)
        CooccurrenceMatrix().compute(self.samples[0], self.session)
        self.session.commit()
        count = self.session.query(VariantCooccurrence).count()
        # size of matrix = n*(n-1)/2 + n (one half of matrix + diagonal)
        self.assertEqual(
            count, (self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE * (self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE - 1) / 2) +
                   self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE)

        self._assert_cooccurrence_matrix(assert_cooccurrence=False)

    def test_all_samples(self):
        count = self.session.query(VariantCooccurrence).count()
        self.assertEqual(count, 0)
        for s in self.samples:
            CooccurrenceMatrix().compute(s, self.session)
            self.session.commit()
        count = self.session.query(VariantCooccurrence).count()
        self.assertLess(
            count,
            self.NUM_SAMPLES * self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE *
            (self.NUM_VARIANT_OBSERVATIONS_PER_SAMPLE - 1) / 2)
        self._assert_cooccurrence_matrix()

    def _assert_cooccurrence_matrix(self, assert_cooccurrence=True):
        variant_cooccurrences = self.session.query(VariantCooccurrence).all()
        found_cooccurrent_variant = False
        for vo in variant_cooccurrences:
            self.assertLessEqual(vo.chromosome_one, vo.chromosome_two)
            if vo.chromosome_one == vo.chromosome_two:
                self.assertLessEqual(vo.position_one, vo.position_two)
            self.assertLessEqual(vo.count, self.NUM_SAMPLES)
            self.assertGreater(vo.count, 0)
            if vo.count > 1:
                found_cooccurrent_variant = True
        if assert_cooccurrence:
            self.assertTrue(found_cooccurrent_variant)

