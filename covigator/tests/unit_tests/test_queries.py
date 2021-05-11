from itertools import combinations
from unittest import TestCase
from faker import Faker
import numpy as np
from covigator.database.database import Database
from covigator.database.model import JobStatus, DataSource, Sample, Gene
from covigator.database.queries import Queries
from covigator.tests.unit_tests.faked_objects import FakeConfiguration
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample, get_mocked_log, get_mocked_variant, \
    get_mocked_variant_cooccurrence, get_mocked_variant_observation


class QueriesTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True, verbose=True, config=FakeConfiguration())
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
            if job.status == JobStatus.FINISHED:
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
            if job.status == JobStatus.FINISHED:
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
        logs = [get_mocked_log(faker=self.faker, source=DataSource.ENA) for _ in range(25)] + \
                [get_mocked_log(faker=self.faker, source=DataSource.GISAID) for _ in range(25)]
        self.session.add_all(logs)
        self.session.commit()
        observed_date = self.queries.get_date_of_last_check(data_source=DataSource.ENA)
        self.assertIsNotNone(observed_date)
        observed_date_gisaid = self.queries.get_date_of_last_check(data_source=DataSource.GISAID)
        self.assertIsNotNone(observed_date_gisaid)

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
        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S", test=True)
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
        variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p):(v1, v2) for v1, v2 in
                              list(combinations(variants, 2))}
        other_variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p): (v1, v2) for v1, v2 in
                              list(combinations(other_variants, 2))}
        combined_variants_to_sample = {"{}-{}".format(v1.hgvs_p, v2.hgvs_p): (v1, v2) for v1, v2 in
                                    list(zip(variants, other_variants))}

        for variant in variants + other_variants:
            cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant, variant))
        for (variant_one, variant_two) in [variants_to_sample.get(k) for k in
                                           np.random.choice(list(variants_to_sample.keys()), 5, replace=False)]:
            cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        for (variant_one, variant_two) in [other_variants_to_sample.get(k) for k in
                                           np.random.choice(list(other_variants_to_sample.keys()), 5, replace=False)]:
            other_cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        for (variant_one, variant_two) in [combined_variants_to_sample.get(k) for k in
                                           np.random.choice(list(combined_variants_to_sample.keys()), 5, replace=False)]:
            other_cooccurrences.append(get_mocked_variant_cooccurrence(self.faker, variant_one, variant_two))
        self.session.add_all(cooccurrences + other_cooccurrences)
        self.session.commit()

        # add some samples to compute the frequency right
        for _ in range(10):
            sample_ena, sample, job = get_mocked_ena_sample(Faker())
            self.session.add(sample_ena)
            self.session.add(sample)
            self.session.add(job)
            self.session.commit()

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S", min_cooccurrence=1, test=True)
        self.assertIsNotNone(data)
        num_unique_variants = len(set([v.variant_id for v in variants]))
        self.assertEqual(data.shape[0], num_unique_variants * num_unique_variants)
        self.assertEqual(data.shape[1], 5)
        self.assertGreater(data[data["count"] > 0].shape[0], len(variants))

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="S", min_cooccurrence=11, test=True)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="X", min_cooccurrence=1, test=True)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_by_gene(gene_name="N", min_cooccurrence=1, test=True)
        self.assertIsNotNone(data)
        self.assertEqual(data.shape[1], 5)
        self.assertGreater(data[data["count"] > 0].shape[0], len(other_variants))

    def test_get_variant_abundance_histogram(self):

        # gets an empty histogram
        histogram = self.queries.get_variant_abundance_histogram()
        self.assertIsNone(histogram)

        # mocks a 100 variants
        num_variants = 100
        variants = [get_mocked_variant(faker=self.faker, chromosome="chr_test") for _ in range(num_variants)]
        self.session.add_all(variants)
        self.session.commit()

        # gets an histogram over 100 variants without variant observations
        histogram = self.queries.get_variant_abundance_histogram()
        self.assertIsNotNone(histogram)
        self.assertGreater(histogram.shape[0], 0)
        self.assertEqual(histogram.shape[1], 3)
        self.assertEqual(histogram.count_unique_variants.sum(), num_variants)
        self.assertEqual(histogram.count_variant_observations.sum(), 0)

        # mock some variant observations
        test_samples = [Sample(id=self.faker.unique.uuid4(), source=DataSource.ENA) for _ in range(5)]
        self.session.add_all(test_samples)
        self.session.commit()
        variant_observations = []
        for v in variants:
            for i in range(self.faker.random_int(min=1, max=5)):
                variant_observations.append(get_mocked_variant_observation(variant=v, sample=test_samples[i]))
        self.session.add_all(variant_observations)
        self.session.commit()

        # gets an histogram over 100 variants
        histogram = self.queries.get_variant_abundance_histogram()
        self.assertIsNotNone(histogram)
        self.assertGreater(histogram.shape[0], 0)
        self.assertEqual(histogram.shape[1], 3)
        self.assertEqual(histogram.count_unique_variants.sum(), num_variants)
        self.assertEqual(histogram.count_variant_observations.sum(), len(variant_observations))

    def test_get_conservation_table(self):
        conservation50 = self.queries.get_conservation_table(bin_size=50)
        self.assertIsNotNone(conservation50)
        self.assertEqual(conservation50.shape[1], 4)
        self.assertGreater(conservation50.shape[0], 0)
        self.assertEqual(conservation50[conservation50.conservation.isna()].shape[0], 0)
        self.assertEqual(conservation50[conservation50.conservation_sarbecovirus.isna()].shape[0], 0)
        self.assertEqual(conservation50[conservation50.conservation_vertebrates.isna()].shape[0], 0)

        conservation10 = self.queries.get_conservation_table(bin_size=10)
        self.assertGreater(conservation10.shape[0], conservation50.shape[0])
        self.assertEqual(conservation10.shape[1], conservation50.shape[1])

        conservation5 = self.queries.get_conservation_table(bin_size=5)
        self.assertGreater(conservation5.shape[0], conservation50.shape[0])
        self.assertGreater(conservation5.shape[0], conservation10.shape[0])
        self.assertEqual(conservation5.shape[1], conservation50.shape[1])

        conservation5_smaller = self.queries.get_conservation_table(bin_size=5, start=10000, end=20000)
        self.assertGreater(conservation5.shape[0], conservation5_smaller.shape[0])
        self.assertEqual(conservation5.shape[1], conservation50.shape[1])

    def test_get_genes_metadata(self):
        genes = self.queries.get_genes_metadata()
        self.assertIsNotNone(genes)
        self.assertGreater(len(genes), 0)
        for g in genes:
            self.assertIsInstance(g, Gene)

    def test_get_genes(self):
        genes = self.queries.get_genes()
        self.assertIsNotNone(genes)
        self.assertGreater(len(genes), 0)
        for g in genes:
            self.assertIsInstance(g, str)

    def test_get_gene(self):
        gene = self.queries.get_gene("S")
        self.assertIsNotNone(gene)
        self.assertIsInstance(gene, Gene)
        gene = self.queries.get_gene("NOEXISTO")
        self.assertIsNone(gene)

    def test_count_jobs_in_queue(self):
        samples = [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.QUEUED) for _ in range(50)] + \
                  [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.FINISHED) for _ in range(50)]
        for sample_ena, sample, job in samples:
            self.session.add_all([sample_ena, sample, job])
        self.session.commit()
        count_jobs_in_queue = self.queries.count_jobs_in_queue(DataSource.ENA)
        self.assertEqual(count_jobs_in_queue, 50)
        count_jobs_in_queue = self.queries.count_jobs_in_queue(DataSource.GISAID)
        self.assertEqual(count_jobs_in_queue, 0)
