from dash_core_components import Markdown, Graph
from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.figures.variants import VariantsFigures
from covigator.database.model import VariantType
from covigator.precomputations.loader import PrecomputationsLoader
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.database.queries import Queries
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, mock_samples


class FiguresTests(AbstractTest):

    def setUp(self) -> None:
        self.sample_figures = SampleFigures(queries=Queries(session=self.session))
        self.variants_figures = VariantsFigures(queries=Queries(session=self.session))

    def test_samples_by_country(self):
        # populates the ENA samples tables
        mock_samples(faker=self.faker, session=self.session, num_samples=100)
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
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        PrecomputationsLoader(session=self.session).load_counts_variants_per_sample()
        figure = self.sample_figures.get_variants_per_sample_plot()
        self.assertIsNotNone(figure)
        self.assertTrue(len(figure) == 2)
        self.assertIsInstance(figure[0], Graph)
        self.assertIsInstance(figure[1], Markdown)

    def test_variants_per_sample_no_data(self):
        figure = self.sample_figures.get_variants_per_sample_plot()
        self.assertIsInstance(figure, Markdown)

    def test_substitutions(self):
        # populates the ENA samples tables
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        PrecomputationsLoader(session=self.session).load_count_substitutions()
        figure = self.sample_figures.get_substitutions_plot(variant_types=[VariantType.SNV.name])
        self.assertIsNotNone(figure)
        self.assertTrue(len(figure) == 2)
        self.assertIsInstance(figure[0], Graph)
        self.assertIsInstance(figure[1], Markdown)

    def test_substitutions_no_data(self):
        figure = self.sample_figures.get_substitutions_plot(variant_types=[VariantType.SNV.name])
        self.assertIsInstance(figure, Markdown)

    def test_needle_plot_no_data(self):
        figure = self.variants_figures.get_variants_plot(gene_name="S", selected_variants=None, bin_size=50)
        self.assertIsInstance(figure, Markdown)

    def test_dn_ds_plot(self):
        figure = self.sample_figures.get_dnds_by_gene_plot()
        self.assertIsInstance(figure, Markdown)

        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        loader = NsSCountsLoader(session=self.session)
        loader.load()
        figure = self.sample_figures.get_dnds_by_gene_plot()
        self.assertIsInstance(figure, list)
        self.assertIsInstance(figure[0], Graph)
        self.assertIsInstance(figure[1], Markdown)
