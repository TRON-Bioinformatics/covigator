from parameterized import parameterized
from dash_core_components import Markdown, Graph
from covigator.dashboard.figures.mutation_stats import MutationStatsFigures
from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.figures.recurrent_mutations import RecurrentMutationsFigures
from covigator.database.model import VariantType, DataSource
from covigator.exceptions import CovigatorQueryException
from covigator.precomputations.loader import PrecomputationsLoader
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.database.queries import Queries
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, mock_samples


class FiguresTests(AbstractTest):

    def setUp(self) -> None:
        self.sample_figures = SampleFigures(queries=Queries(session=self.session))
        self.mutation_stats_figures = MutationStatsFigures(queries=Queries(session=self.session))
        self.variants_figures = RecurrentMutationsFigures(queries=Queries(session=self.session))

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_samples_by_country(self, source):
        # populates the ENA samples tables
        mock_samples(faker=self.faker, session=self.session, num_samples=100)
        if source is not None:
            figure = self.sample_figures.get_accumulated_samples_by_country_plot(
                min_samples=0, data_source=source)
            self.assertIsNotNone(figure)
            self.assertTrue(len(figure) == 2)
            self.assertIsInstance(figure[0], Graph)
            self.assertIsInstance(figure[1], Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.sample_figures.get_accumulated_samples_by_country_plot,
                data_source=source
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_samples_by_country_no_data(self, source):
        if source is not None:
            figure = self.sample_figures.get_accumulated_samples_by_country_plot(data_source=source)
            self.assertIsInstance(figure, Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.sample_figures.get_accumulated_samples_by_country_plot,
                data_source=source
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_variants_per_sample(self, source):
        # populates the ENA samples tables
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        PrecomputationsLoader(session=self.session).load_counts_variants_per_sample()
        if source is not None:
            figure = self.mutation_stats_figures.get_variants_per_sample_plot(data_source=source)
            self.assertIsNotNone(figure)
            self.assertTrue(len(figure) == 2)
            self.assertIsInstance(figure[0], Graph)
            self.assertIsInstance(figure[1], Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.mutation_stats_figures.get_variants_per_sample_plot,
                data_source=source
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_variants_per_sample_no_data(self, source):
        if source is not None:
            figure = self.mutation_stats_figures.get_variants_per_sample_plot(data_source=source)
            self.assertIsInstance(figure, Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.mutation_stats_figures.get_variants_per_sample_plot,
                data_source=source
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_substitutions(self, source):
        # populates the ENA samples tables
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        PrecomputationsLoader(session=self.session).load_count_substitutions()

        if source is not None:
            figure = self.mutation_stats_figures.get_substitutions_plot(
                variant_types=[VariantType.SNV.name], data_source=source)
            self.assertIsNotNone(figure)
            self.assertTrue(len(figure) == 2)
            self.assertIsInstance(figure[0], Graph)
            self.assertIsInstance(figure[1], Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.mutation_stats_figures.get_substitutions_plot,
                data_source=source, variant_types=[VariantType.SNV.name]
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_substitutions_no_data(self, source):
        if source is not None:
            figure = self.mutation_stats_figures.get_substitutions_plot(
                variant_types=[VariantType.SNV.name], data_source=source)
            self.assertIsInstance(figure, Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.mutation_stats_figures.get_substitutions_plot,
                data_source=source, variant_types=[VariantType.SNV.name]
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_needle_plot(self, source):
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        if source is not None:
            figure = self.variants_figures.get_variants_plot(
                gene_name="S", domain_name=None, selected_variants=None, bin_size=50, source=source)
            self.assertIsInstance(figure, Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.variants_figures.get_variants_plot,
                source=source, gene_name="S", bin_size=50, domain_name=None, selected_variants=None
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_needle_plot_no_data(self, source):
        if source is not None:
            figure = self.variants_figures.get_variants_plot(
                gene_name="S", domain_name=None, selected_variants=None, bin_size=50, source=source)
            self.assertIsInstance(figure, Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.variants_figures.get_variants_plot,
                source=source, gene_name="S", bin_size=50, domain_name=None, selected_variants=None
            )

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name, (None,)])
    def test_dn_ds_plot(self, source):
        if source is not None:
            figure = self.sample_figures.get_dnds_by_gene_plot(data_source=source)
            self.assertIsInstance(figure, Markdown)

        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        loader = NsSCountsLoader(session=self.session)
        loader.load()

        if source is not None:
            figure = self.sample_figures.get_dnds_by_gene_plot(data_source=source)
            self.assertIsInstance(figure, list)
            self.assertIsInstance(figure[0], Graph)
            self.assertIsInstance(figure[1], Markdown)
        else:
            self.assertRaises(
                CovigatorQueryException,
                self.sample_figures.get_dnds_by_gene_plot,
                data_source=source
            )
