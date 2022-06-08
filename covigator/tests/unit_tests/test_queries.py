import unittest

from parameterized import parameterized

from covigator import SYNONYMOUS_VARIANT
from covigator.precomputations.loader import PrecomputationsLoader
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.database.model import JobStatus, DataSource, Gene, RegionType, Domain
from covigator.database.queries import Queries
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import get_mocked_variant, \
    get_mocked_variant_observation, mock_samples, mock_cooccurrence_matrix, mock_samples_and_variants, MOCKED_DOMAINS, \
    get_mocked_sample


class QueriesTests(AbstractTest):

    def setUp(self) -> None:
        self.queries = Queries(session=self.session)

    def test_get_date_of_first_ena_sample(self):
        first_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = mock_samples(faker=self.faker, session=self.session, num_samples=50, source=DataSource.ENA.name) + \
                  mock_samples(faker=self.faker, session=self.session, job_status=JobStatus.FAILED_PROCESSING, num_samples=50,
                               source=DataSource.ENA.name)
        for sample_ena in samples:
            if sample_ena.status == JobStatus.FINISHED:
                if first_sample_date is None:
                    first_sample_date = sample_ena.collection_date
                if sample_ena.collection_date < first_sample_date:
                    first_sample_date = sample_ena.collection_date
        observed_date = self.queries.get_date_of_first_sample()
        self.assertEqual(observed_date, first_sample_date)

    def test_get_date_of_first_ena_sample_empty(self):
        observed_date = self.queries.get_date_of_first_sample()
        self.assertIsNone(observed_date)

    def test_get_date_of_most_recent_ena_sample(self):
        most_recent_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = mock_samples(faker=self.faker, session=self.session, num_samples=50, source=DataSource.ENA.name) + \
                  mock_samples(faker=self.faker, session=self.session, job_status=JobStatus.FAILED_PROCESSING, num_samples=50,
                               source=DataSource.ENA.name)
        for sample_ena in samples:
            if sample_ena.status == JobStatus.FINISHED:
                if most_recent_sample_date is None:
                    most_recent_sample_date = sample_ena.collection_date
                if sample_ena.collection_date > most_recent_sample_date:
                    most_recent_sample_date = sample_ena.collection_date
        observed_date = self.queries.get_date_of_most_recent_sample()
        self.assertEqual(observed_date, most_recent_sample_date)

    def test_get_date_of_most_recent_ena_sample_empty(self):
        observed_date = self.queries.get_date_of_most_recent_sample()
        self.assertIsNone(observed_date)

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name])
    def test_get_cooccurrence_matrix_by_gene_no_data(self, source):
        PrecomputationsLoader(session=self.session).load_table_counts()
        data = self.queries.get_sparse_cooccurrence_matrix(gene_name="S", domain=None, source=source)
        self.assertEqual(data.shape, (0, 10))

    # TODO: make this test stable
    @unittest.skip
    def test_get_cooccurrence_matrix_by_gene(self):
        other_variants, variants = mock_cooccurrence_matrix(faker=self.faker, session=self.session)
        PrecomputationsLoader(session=self.session).load_table_counts()

        data = self.queries.get_variants_cooccurrence_matrix(gene_name="S", min_cooccurrence=1, test=True)
        self.assertIsNotNone(data)
        num_unique_variants = len(set([v.variant_id for v in variants if v.annotation != SYNONYMOUS_VARIANT]))
        self.assertEqual(data.shape[0], num_unique_variants * num_unique_variants)
        self.assertEqual(data.shape[1], 7)
        self.assertGreaterEqual(data[data["count"] > 0].shape[0], len(variants))
        self.assertGreaterEqual(data[data["frequency"] > 0].shape[0], len(variants))
        self.assertGreaterEqual(data[data["jaccard"] > 0].shape[0], len(variants))

        data = self.queries.get_variants_cooccurrence_matrix(gene_name="S", min_cooccurrence=11, test=True)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_matrix(gene_name="X", min_cooccurrence=1, test=True)
        self.assertIsNone(data)

        data = self.queries.get_variants_cooccurrence_matrix(gene_name="N", min_cooccurrence=1, test=True)
        self.assertIsNotNone(data)
        self.assertEqual(data.shape[1], 7)
        self.assertGreaterEqual(data[data["count"] > 0].shape[0], len(other_variants))
        self.assertGreaterEqual(data[data["frequency"] > 0].shape[0], len(variants))
        self.assertGreaterEqual(data[data["jaccard"] > 0].shape[0], len(variants))

    @unittest.skip
    def test_get_mds(self):
        mock_cooccurrence_matrix(faker=self.faker, session=self.session)

        mds_fit, mds_coords = self.queries.get_mds(gene_name="S")
        self.assertIsNotNone(mds_fit)

    @parameterized.expand([(DataSource.ENA, ), (DataSource.GISAID, )])
    def test_get_variant_abundance_histogram(self, source):

        # gets an empty histogram
        histogram = self.queries.get_variant_abundance_histogram(cache=False, source=source.name)
        self.assertIsNone(histogram)

        # mocks a 100 variants
        num_variants = 100
        variants = [get_mocked_variant(faker=self.faker, chromosome="chr_test", source=source.name) for _ in range(num_variants)]
        self.session.add_all(variants)
        self.session.commit()

        # gets an histogram over 100 variants without variant observations
        histogram = self.queries.get_variant_abundance_histogram(cache=False, source=source.name)
        self.assertIsNotNone(histogram)
        self.assertGreater(histogram.shape[0], 0)
        self.assertEqual(histogram.shape[1], 3)
        self.assertEqual(histogram.count_unique_variants.sum(), num_variants)
        self.assertEqual(histogram.count_variant_observations.sum(), 0)

        # mock some variant observations
        test_samples = [get_mocked_sample(faker=self.faker, source=source) for _ in range(5)]
        self.session.add_all(test_samples)
        self.session.commit()
        variant_observations = []
        for v in variants:
            for i in range(self.faker.random_int(min=1, max=5)):
                variant_observations.append(get_mocked_variant_observation(variant=v, sample=test_samples[i]))
        self.session.add_all(variant_observations)
        self.session.commit()

        # gets an histogram over 100 variants
        histogram = self.queries.get_variant_abundance_histogram(cache=False, source=source.name)
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
        genes = self.queries.get_genes()
        self.assertIsNotNone(genes)
        self.assertGreater(len(genes), 0)
        for g in genes:
            self.assertIsInstance(g, Gene)

    def test_get_genes(self):
        genes = self.queries.get_genes()
        self.assertIsNotNone(genes)
        self.assertGreater(len(genes), 0)
        for g in genes:
            self.assertIsInstance(g, Gene)
            self.assertIsNotNone(g.name)

    def test_get_gene(self):
        gene = self.queries.get_gene("S")
        self.assertIsNotNone(gene)
        self.assertIsInstance(gene, Gene)
        gene = self.queries.get_gene("NOEXISTO")
        self.assertIsNone(gene)

    def test_get_domains(self):
        domains = self.queries.get_domains()
        self.assertIsNotNone(domains)
        self.assertGreater(len(domains), 0)
        for d in domains:
            self.assertIsInstance(d, Domain)
            self.assertIsNotNone(d.name)
            self.assertIsNotNone(d.gene_name)

    def test_get_domains_by_gene(self):
        domains = self.queries.get_domains_by_gene("S")
        self.assertIsNotNone(domains)
        self.assertGreater(len(domains), 0)
        for d in domains:
            self.assertIsInstance(d, Domain)
            self.assertIsNotNone(d.name)
            self.assertEqual(d.gene_name, "S")

    def test_count_jobs_in_queue(self):
        mock_samples(faker=self.faker, session=self.session, job_status=JobStatus.QUEUED, num_samples=50)
        mock_samples(faker=self.faker, session=self.session, job_status=JobStatus.FINISHED, num_samples=50)

        count_jobs_in_queue_ena = self.queries.count_jobs_in_queue(DataSource.ENA)
        self.assertGreater(count_jobs_in_queue_ena, 0)
        count_jobs_in_queue_gisaid = self.queries.count_jobs_in_queue(DataSource.GISAID)
        self.assertGreater(count_jobs_in_queue_gisaid, 0)

        self.assertEqual(count_jobs_in_queue_ena + count_jobs_in_queue_gisaid, 50)

    def test_get_dnds_table(self):
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        nsSCountsLoader = NsSCountsLoader(session=self.session)
        nsSCountsLoader.load()

        data = self.queries.get_dnds_table(source=DataSource.ENA.name)
        self._assert_dnds_table(data)
        self.assertEqual(data[data.source != DataSource.ENA].shape[0], 0)      # no entries to other source

        data = self.queries.get_dnds_table(source=DataSource.GISAID.name)
        self._assert_dnds_table(data)
        self.assertEqual(data[data.source != DataSource.GISAID].shape[0], 0)  # no entries to other source

        countries = list(data.country.unique())[0:2]
        data = self.queries.get_dnds_table(source=DataSource.GISAID.name, countries=countries)
        self._assert_dnds_table(data)
        self.assertEqual(data[~data.country.isin(countries)].shape[0], 0)  # no entries to other country

        genes = ["S"]
        data = self.queries.get_dnds_table(source=DataSource.ENA.name, genes=genes)
        self._assert_dnds_table(data, has_domains=True)
        self.assertEqual(data[~data.region_name.isin(genes + MOCKED_DOMAINS)].shape[0], 0)  # no entries to other country

    def _assert_dnds_table(self, data, has_domains=False):
        self.assertIsNotNone(data)
        self.assertGreater(data.shape[0], 0)
        self.assertEqual(data.shape[1], 8)
        self.assertGreater(data[data.region_type == RegionType.GENE].shape[0], 0)
        if has_domains:
            self.assertGreater(data[data.region_type == RegionType.DOMAIN].shape[0], 0)
        self.assertEqual(data[(data.region_type != RegionType.GENE) & (data.region_type != RegionType.DOMAIN)].shape[0], 0)  # all entries to a gene or domain
        self.assertEqual(data[data.country.isna()].shape[0], 0)  # no empty countries
        self.assertEqual(data[data.month.isna()].shape[0], 0)  # no empty months
        self.assertEqual(data[data.region_name.isna()].shape[0], 0)  # no empty gene
