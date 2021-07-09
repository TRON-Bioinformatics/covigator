from math import sqrt
from typing import List
import numpy as np
import pandas as pd
import plotly
from logzero import logger

from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, MARGIN, TEMPLATE
import plotly.express as px
import plotly.graph_objects as go
import dash_core_components as dcc

from covigator.database.model import VariantType, DataSource


INDEL_TYPE_COLOR_MAP = {
    "INSERTION_INFRAME": "#fc8d62",
    "DELETION_INFRAME": "#66c2a5",
    "INSERTION_FRAMESHIFT": "#fdc4ad",
    "DELETION_FRAMESHIFT": "#9dd8c5",
}


class SampleFigures(Figures):

    def get_annotations_plot(self, data_source: str = None, genes: List[str] = None):
        data = self.queries.get_annotations(data_source=data_source, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            fig = px.bar(data, y="count", x="annotation", log_y=True, color='count',
                         color_continuous_scale=plotly.colors.sequential.Brwnyl)
            #.update_yaxes(categoryorder="total descending")
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                yaxis={'title': "num. variants (log)"},  # , 'autorange': 'reversed'},
                xaxis={'title': None},
            )
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Most common mutation effects**
                
                *Ratio of non synonymous to synonymous SNVs (dN/dS): {dnds}*
                """.format(dnds=round(data[data.annotation == "missense_variant"]["count"].sum() /
                                data[data.annotation == "synonymous_variant"]["count"].sum(), 3)))
            ]
        return graph

    def get_indels_lengths_plot(self, data_source: str = None, genes: List[str] = None):
        data = self.queries.get_indel_lengths(data_source=data_source, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            fig = px.bar(
                data, y="count", x="length", color="variant_type", color_discrete_map=INDEL_TYPE_COLOR_MAP)
                #.update_yaxes(categoryorder="total descending")
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                yaxis={'title': "num. samples"},  #, 'autorange': 'reversed'},
                xaxis={'title': None},
            )
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Indel length distribution**
                
                *Insertion to deletion ratio: {indel_ratio}*
                """.format(
                    indel_ratio=round(data[
                        data.variant_type.isin(["INSERTION_INFRAME", "INSERTION_FRAMESHIFT"])]["count"].sum() /
                        data[data.variant_type.isin(["DELETION_INFRAME", "DELETION_FRAMESHIFT"])]["count"].sum(), 3)))
            ]
        return graph

    def get_substitutions_plot(self, variant_types: List[str], data_source: str = None, genes: List[str] = None):
        data = self.queries.get_substitutions(data_source=data_source, genes=genes, variant_types=variant_types)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            fig = px.bar(
                data, y="substitution", x="count", color="variant_type", text="rate",
                color_discrete_map=VARIANT_TYPE_COLOR_MAP)\
                .update_yaxes(categoryorder="total descending")
            fig.update_traces(textposition='outside')
            fig.update_layout(
                margin=go.layout.Margin(l=0, r=40, b=0, t=30),  # we need some extra space for some labels overflowing
                template=TEMPLATE,
                showlegend=False,
                yaxis={'title': None, 'autorange': 'reversed'},
                xaxis={'title': "num. samples"},
                uniformtext_mode='show',
                uniformtext_minsize=8
            )
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Top mutations**
                
                *Only mutations occurring at least in 10 samples are represented*
                """)
            ]
        return graph

    def get_variants_per_sample_plot(
            self, data_source: str = None, genes: List[str] = None, variant_types: List[str] = None):

        data = self.queries.get_variants_per_sample(data_source=data_source, genes=genes, variant_types=variant_types)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            counts = np.repeat(data.number_mutations, data["count"])
            fig = px.violin(y=counts, orientation='v', box=True, color_discrete_sequence=["grey"])
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                yaxis={'title': "num. mutations"},
                xaxis={'title': None},
            )

            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                        **Mutations per sample**

                        *Median: {median} (IQR: {iqr})*
                        """.format(median=round(np.median(counts), 3),
                                   iqr=round(np.percentile(counts, 75) - np.percentile(counts, 25), 3)))
            ]
        return graph

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
                          color_discrete_sequence=px.colors.qualitative.Light24)
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

                *Top {} countries: {}. Countries with < {} samples are excluded*
                """.format(len(top_countries_tooltip),
                           ", ".join(top_countries_tooltip),
                           min_samples))
            ]
        return graph
