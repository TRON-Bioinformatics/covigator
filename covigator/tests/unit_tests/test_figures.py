from unittest import TestCase

from dash_core_components import Markdown, Graph
from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.figures.variants import VariantsFigures
from covigator.database.database import Database
from faker import Faker

from covigator.database.precomputed import Precomputer
from covigator.database.queries import Queries
from covigator.tests.unit_tests.faked_objects import FakeConfiguration
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample, get_mocked_variant, get_mocked_variant_observation


class FiguresTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True, config=FakeConfiguration())
        self.session = self.database.get_database_session()
        self.sample_figures = SampleFigures(queries=Queries(session=self.session))
        self.variants_figures = VariantsFigures(queries=Queries(session=self.session))
        self.faker = Faker()

    def test_samples_by_country(self):
        # populates the ENA samples tables
        for _ in range(100):
            sample_ena, sample, job = get_mocked_ena_sample(faker=self.faker)
            self.session.add_all([sample_ena, sample, job])
        self.session.commit()
        figure = self.sample_figures.get_accumulated_samples_by_country_plot(min_samples=0)
        self.assertIsNotNone(figure)
        self.assertTrue(len(figure) == 2)
        self.assertIsInstance(figure[0], Graph)
        self.assertIsInstance(figure[1], Markdown)

    def test_samples_by_country_no_data(self):
        figure = self.sample_figures.get_accumulated_samples_by_country_plot()
        self.assertIsInstance(figure, Markdown)

    def test_variants_per_sample(self):
        # populates the ENA samples tables
        existing_variants = set()
        for _ in range(100):
            sample_ena, sample, job = get_mocked_ena_sample(faker=self.faker)
            self.session.add_all([sample_ena, sample, job])
            variants = [get_mocked_variant(faker=self.faker) for _ in range(10)]
            self.session.add_all(list(filter(lambda x: x in existing_variants, variants)))
            existing_variants.update([v.variant_id for v in variants])
            self.session.commit()
            variants_observations = [get_mocked_variant_observation(faker=self.faker, variant=v, sample=sample)
                                     for v in variants]
            self.session.add_all(variants_observations)
        self.session.commit()
        Precomputer(session=self.session).load_counts_variants_per_sample()
        figure = self.sample_figures.get_variants_per_sample_plot()
        self.assertIsNotNone(figure)
        self.assertTrue(len(figure) == 2)
        self.assertIsInstance(figure[0], Graph)
        self.assertIsInstance(figure[1], Markdown)

    def test_variants_per_sample_no_data(self):
        figure = self.sample_figures.get_variants_per_sample_plot()
        self.assertIsInstance(figure, Markdown)

    def test_needle_plot_no_data(self):
        figure = self.variants_figures.get_variants_plot(gene_name="S", selected_variants=None, bin_size=50)
        self.assertIsInstance(figure, Markdown)
