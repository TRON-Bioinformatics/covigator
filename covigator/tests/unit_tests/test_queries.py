import random
from itertools import combinations
from unittest import TestCase, skip
from faker import Faker
from covigator.database.database import Database
from covigator.database.model import JobStatus, DataSource
from covigator.database.queries import Queries
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample, get_mocked_log, get_mocked_variant, \
    get_mocked_variant_cooccurrence


class QueriesTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True, verbose=True)
        self.session = self.database.get_database_session()
        self.queries = Queries(session=self.session)
        self.faker = Faker()

    def test_get_date_of_first_ena_sample(self):
        first_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = [get_mocked_ena_sample(faker=self.faker) for _ in range(50)] + \
                  [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.FAILED_LOAD) for _ in range(50)]
        for sample_ena, sample, job in samples:
            self.session.add_all([sample_ena, sample, job])
            if job.status == JobStatus.COOCCURRENCE:
                if first_sample_date is None:
                    first_sample_date = sample_ena.first_created
                if sample_ena.first_created < first_sample_date:
                    first_sample_date = sample_ena.first_created
        self.session.commit()
        observed_date = self.queries.get_date_of_first_ena_sample()
        self.assertEqual(observed_date, first_sample_date.date())

    def test_get_date_of_first_ena_sample_empty(self):
        observed_date = self.queries.get_date_of_first_ena_sample()
        self.assertIsNone(observed_date)

    def test_get_date_of_most_recent_ena_sample(self):
        most_recent_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = [get_mocked_ena_sample(faker=self.faker) for _ in range(50)] + \
                  [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.FAILED_LOAD) for _ in range(50)]
        for sample_ena, sample, job in samples:
            self.session.add_all([sample_ena, sample, job])
            if job.status == JobStatus.COOCCURRENCE:
                if most_recent_sample_date is None:
                    most_recent_sample_date = sample_ena.first_created
                if sample_ena.first_created > most_recent_sample_date:
                    most_recent_sample_date = sample_ena.first_created
        self.session.commit()
        observed_date = self.queries.get_date_of_most_recent_ena_sample()
        self.assertEqual(observed_date, most_recent_sample_date.date())

    def test_get_date_of_most_recent_ena_sample_empty(self):
        observed_date = self.queries.get_date_of_most_recent_ena_sample()
        self.assertIsNone(observed_date)

    def test_get_date_of_last_ena_check(self):
        logs = [get_mocked_log(faker=self.faker) for _ in range(50)]
        self.session.add_all(logs)
        self.session.commit()
        observed_date = self.queries.get_date_of_last_check(data_source=DataSource.ENA)
        self.assertIsNotNone(observed_date)

    def test_get_date_of_last_ena_check_empty(self):
        observed_date = self.queries.get_date_of_last_check(data_source=DataSource.ENA)
        self.assertIsNone(observed_date)
        observed_date = self.queries.get_date_of_last_check(data_source=DataSource.GISAID)
        self.assertIsNone(observed_date)

    def test_get_date_of_last_ena_update(self):
        # NOTE: the implementation that works in Postgres does not work in SQLite!!
        logs = [get_mocked_log(faker=self.faker, source=DataSource.ENA) for _ in range(50)]
        self.session.add_all(logs)
        self.session.commit()
        observed_date = self.queries.get_date_of_last_update(data_source=DataSource.ENA)
        self.assertIsNotNone(observed_date)

    def test_get_date_of_last_ena_update_empty(self):
        observed_date = self.queries.get_date_of_last_update(data_source=DataSource.ENA)
        self.assertIsNone(observed_date)
        observed_date = self.queries.get_date_of_last_update(data_source=DataSource.GISAID)
        self.assertIsNone(observed_date)

    def test_get_cooccurrence_matrix_by_gene_no_data(self):
        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S")
        self.assertIsNone(data)

    def test_get_cooccurrence_matrix_by_gene(self):
        # add some variants belonging to two genes
        chromosome = "fixed_chromosome"
        gene_name = "S"
        variants = [get_mocked_variant(faker=self.faker, chromosome=chromosome, gene_name=gene_name) for _ in range(5)]
        other_gene_name = "N"
        other_variants = [get_mocked_variant(faker=self.faker, chromosome=chromosome, gene_name=other_gene_name) for _ in range(5)]
        self.session.add_all(variants + other_variants)
        self.session.commit()

        # adds some cooccurrences
        cooccurrences = []
        other_cooccurrences = []
        for (variant_one, variant_two) in random.sample(list(combinations(variants, 2)), 5):
            cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        for (variant_one, variant_two) in random.sample(list(combinations(other_variants, 2)), 5):
            other_cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        for (variant_one, variant_two) in random.sample(list(zip(variants, other_variants)), 5):
            other_cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        self.session.add_all(cooccurrences + other_cooccurrences)
        self.session.commit()

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S", min_cooccurrence=1)
        self.assertIsNotNone(data)
        num_unique_variants = len(set([c.position_one for c in cooccurrences] + [c.position_two for c in cooccurrences]))
        self.assertEqual(data.shape[0], num_unique_variants * (num_unique_variants - 1) + num_unique_variants)
        self.assertEqual(data.shape[1], 7)
        self.assertEqual(data[data["count"] > 0].shape[0], 5)

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S", min_cooccurrence=11)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="X", min_cooccurrence=1)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="N", min_cooccurrence=1)
        self.assertIsNotNone(data)
        self.assertEqual(data.shape[1], 7)
        self.assertEqual(data[data["count"] > 0].shape[0], 5)
