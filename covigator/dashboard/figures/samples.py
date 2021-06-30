from math import sqrt
from typing import List

from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, MARGIN, TEMPLATE
import plotly.express as px
import plotly.graph_objects as go
import dash_core_components as dcc

from covigator.database.model import VariantType, DataSource


class SampleFigures(Figures):

    def get_variants_per_sample_plot(self, data_source: str = None, genes: List[str] = None):
        data = self.queries.get_variants_per_sample(data_source=data_source, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            fig = px.bar(data, x="number_mutations", y="count", color="variant_type",
                         color_discrete_map={
                             VariantType.SNV.name: "#8da0cb",
                             VariantType.INSERTION.name: "#fc8d62",
                             VariantType.DELETION.name: "#66c2a5",
                         })
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                xaxis={'title': "num. mutations"},
                yaxis={'title': "num. samples"},
            )
            average, std = self._get_avg_and_std(data[data.variant_type == VariantType.SNV.name])
            average_insertions, std_insertions = self._get_avg_and_std(
                data[data.variant_type == VariantType.INSERTION.name])
            average_deletions, std_deletions = self._get_avg_and_std(
                data[data.variant_type == VariantType.DELETION.name])
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                        **Overall mutations per sample.**

                        *Average SNVs: {average} (STD: {std}). *
                        *Average insertions: {average_insertions} (STD: {std_insertions}). *
                        *Average deletions: {average_deletions} (STD: {std_deletions}).*
                        """.format(average=average, std=std, average_insertions=average_insertions,
                                   std_insertions=std_insertions, average_deletions=average_deletions,
                                   std_deletions=std_deletions))
            ]
        return graph

    def _get_avg_and_std(self, data):
        average = 0.0
        std = 0.0
        if data is not None and data.shape[0] > 0:
            average = round((data["number_mutations"] * data["count"]).sum() / data["count"].sum(), 3)
            std = round(sqrt((data["number_mutations"].transform(
                lambda x: (x - average) ** 2) * data["count"]).sum() / data["count"].sum()), 3)
        return average, std

    def get_accumulated_samples_by_country_plot(self, data_source: DataSource = None, countries=None, min_samples=1000):
        data = self.queries.get_accumulated_samples_by_country(
            data_source=data_source, countries=countries, min_samples=min_samples)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            countries = list(data.sort_values("cumsum", ascending=False).country.unique())
            fig = px.area(data, x="date", y="cumsum", color="country",
                          category_orders={
                              "country": countries[::-1]},
                          labels={"cumsum": "num. samples", "count": "increment"},
                          hover_data=["count"],
                          color_discrete_sequence=px.colors.qualitative.Dark24)
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'reversed', 'title': None},
                xaxis={'title': None},
            )

            top_countries = countries[0: min(5, len(countries))]
            top_countries_and_cumsum = data[data.country.isin(top_countries)][["country", "cumsum"]] \
                .groupby("country").max().reset_index().sort_values("cumsum", ascending=False)
            top_countries_tooltip = list(
                top_countries_and_cumsum.apply(lambda x: "{} ({})".format(x[0], int(x[1])), axis=1))

            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Accumulated samples by country**

                *Top {} countries: {}*
                """.format(len(top_countries_tooltip), ", ".join(top_countries_tooltip)))
            ]
        return graph
