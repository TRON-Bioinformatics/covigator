from unittest import TestCase

from dash_core_components import Markdown

from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.figures.variants import VariantsFigures
from covigator.database.database import Database
from faker import Faker
from covigator.database.queries import Queries
from covigator.tests.unit_tests.faked_objects import FakeConfiguration
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample


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
        figure = self.sample_figures.get_accumulated_samples_by_country_plot()
        self.assertIsNotNone(figure)

    def test_samples_by_country_no_data(self):
        figure = self.sample_figures.get_accumulated_samples_by_country_plot()
        self.assertIsInstance(figure, Markdown)

    def test_needle_plot_no_data(self):
        figure = self.variants_figures.get_variants_plot(gene_name="S", selected_variants=None, bin_size=50)
        self.assertIsInstance(figure, Markdown)
