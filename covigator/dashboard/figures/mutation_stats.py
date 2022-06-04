from typing import List
import numpy as np
import plotly
from logzero import logger

from covigator import MISSENSE_VARIANT, SYNONYMOUS_VARIANT
from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, MARGIN, TEMPLATE
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc


INDEL_TYPE_COLOR_MAP = {
    "INSERTION_INFRAME": "#ee6002",
    "DELETION_INFRAME": "#09af00",
    "INSERTION_FRAMESHIFT": "#ffddb0",
    "DELETION_FRAMESHIFT": "#defabb",
}


class MutationStatsFigures(Figures):

    def get_annotations_plot(self, data_source: str = None, genes: List[str] = None):
        logger.debug("Getting data on annotations...")
        data = self.queries.get_annotations(data_source=data_source, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Prepare plot on annotations...")
            fig = px.bar(data, y="count", x="annotation", log_y=True, color='count',
                         color_continuous_scale=plotly.colors.sequential.Brwnyl)
            #.update_yaxes(categoryorder="total descending")
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                yaxis={'title': "num. mutations (log)"},  # , 'autorange': 'reversed'},
                xaxis={'title': None},
            )
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Most common mutation effects**
                
                *Ratio of non synonymous to synonymous SNVs (N/S): {dnds}*
                """.format(dnds=round(data[data.annotation == MISSENSE_VARIANT]["count"].sum() /
                                data[data.annotation == SYNONYMOUS_VARIANT]["count"].sum(), 3)))
            ]
        return graph

    def get_indels_lengths_plot(self, data_source: str = None, genes: List[str] = None):
        logger.debug("Getting data on indels lengths...")
        data = self.queries.get_indel_lengths(data_source=data_source, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Prepare plot on indels lengths...")
            fig = px.bar(
                data, y="count", x="length", color="variant_type", color_discrete_map=INDEL_TYPE_COLOR_MAP)
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'normal', 'title': None},
                yaxis={'title': "num. samples"},
                xaxis={'title': "indel length (bp)"}
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

    def get_substitutions_plot(self, data_source: str, variant_types: List[str], genes: List[str] = None):
        logger.debug("Getting data on substitutions plot...")
        data = self.queries.get_substitutions(data_source=data_source, genes=genes, variant_types=variant_types)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Preparing plot on samples by country...")
            fig = px.bar(
                data, y="substitution", x="count", color="variant_type", text="rate",
                color_discrete_map=VARIANT_TYPE_COLOR_MAP)\
                .update_yaxes(categoryorder="total descending")
            fig.update_traces(textposition='outside')
            fig.update_layout(
                margin=go.layout.Margin(l=0, r=40, b=0, t=30),  # we need some extra space for some labels overflowing
                template=TEMPLATE,
                legend={'title': None, 'yanchor': "bottom", 'y': 0.01, 'xanchor': "right", 'x': 0.99},
                showlegend=True,
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
            self, data_source: str, genes: List[str] = None, variant_types: List[str] = None):

        logger.debug("Getting data on variants per sample...")
        data = self.queries.get_variants_per_sample(data_source=data_source, genes=genes, variant_types=variant_types)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Preparing plot on samples by country...")
            counts = np.repeat(data.number_mutations, data["count"])
            median = round(np.median(counts), 3)
            third_quartile = np.percentile(counts, 75)
            first_quartile = np.percentile(counts, 25)
            extreme_threshold = median + (3 * (third_quartile - first_quartile))
            filtered_data = data[data.number_mutations <= extreme_threshold]
            if filtered_data.shape[0] > 0:
                fig = px.bar(filtered_data,
                             x="number_mutations",
                             y='count',
                             color="variant_type",
                             color_discrete_map=VARIANT_TYPE_COLOR_MAP)
                fig.add_vline(x=median, line_width=2, line_dash="dash", line_color='grey',
                              annotation_text="median", annotation_position="top right")
                fig.add_vrect(x0=first_quartile, x1=third_quartile,
                              fillcolor="grey", opacity=0.25, line_width=0)
                fig.update_layout(
                    margin=MARGIN,
                    template=TEMPLATE,
                    legend={'traceorder': 'normal', 'title': None},
                    yaxis={'title': "num. samples"},
                    xaxis={'title': "num. mutations"},
                )

                graph = [
                    dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                    dcc.Markdown("""
                            **Mutations per sample**
    
                            *Median: {median} (IQR: {iqr})*
                            """.format(median=median, iqr=third_quartile - first_quartile))
                ]
        return graph
