from typing import List
import numpy as np
import pandas as pd
from logzero import logger
from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, MARGIN, TEMPLATE
import plotly.express as px
from dash import dcc

from covigator.database.model import DataSource


class SampleFigures(Figures):

    def get_accumulated_samples_by_country_plot(self, data_source: str, countries=None, min_samples=1000,
                                                lineages=None):
        logger.debug("Getting data on samples by country...")
        data = self.queries.get_accumulated_samples_by_country(
            data_source=data_source, countries=countries, min_samples=min_samples, lineages=lineages)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Prepare plot on samples by country...")
            countries = list(data.sort_values("cumsum", ascending=False).country.unique())
            fig = px.area(data, x="date", y="cumsum", color="country",
                          category_orders={
                              "country": countries[::-1]},
                          labels={"cumsum": "num. samples", "count": "increment"},
                          hover_data=["count"],
                          color_discrete_sequence=px.colors.qualitative.Vivid)
            fig.update_traces(line=dict(width=0.5))
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

    def _calculate_dn_ds(self, NS, S, ns, s):
        pn = float(ns) / float(NS)
        ps = float(s) / float(S)
        dn = np.log(1 + pn) + 1
        ds = np.log(1 + ps) + 1
        return dn / ds

    def get_dnds_by_gene_plot(
            self, data_source: str, countries: List[str] =None, genes: List[str] = None):

        logger.debug("Getting data on dN/dS...")
        data = self.queries.get_dnds_table(
            source=data_source, countries=countries, genes=genes)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Prepare plot on dN/dS...")
            genes = pd.concat([
                self.queries.get_genes_df(),
                self.queries.get_domains_df()
            ])
            # prepares the data and calculates the dN/dS
            data_to_plot = data.groupby(["month", "region_name"]).sum().reset_index().sort_values("month")
            data_to_plot = pd.merge(left=genes, right=data_to_plot, left_on="name", right_on="region_name")
            data_to_plot["dn_ds"] = data_to_plot[["ns", "s", "fraction_non_synonymous", "fraction_synonymous"]].apply(
                lambda x: self._calculate_dn_ds(ns=x[0], s=x[1], NS=x[2], S=x[3]), axis=1)

            fig = px.line(data_to_plot, x='month', y='dn_ds', color='region_name',
                          symbol='region_name', line_dash='region_name', line_dash_sequence=['dash'],
                          labels={"dn_ds": "dN/dS", "region_name": "gene"},
                          hover_data=["region_name", "dn_ds"],
                          color_discrete_sequence=px.colors.qualitative.Vivid)
            fig.update_traces(line=dict(width=0.5), marker=dict(size=10))
            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'title': None},
                xaxis={'title': None},
            )

            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **dN/dS by gene**
                """)
            ]
        return graph
