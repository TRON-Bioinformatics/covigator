from typing import List

import colorlover
from dash import dash_table
import pandas as pd
from logzero import logger
from plotly.subplots import make_subplots
from sqlalchemy import and_

from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, MARGIN, TEMPLATE, STYLES_STRIPPED, STYLE_HEADER, \
    STYLE_CELL
import plotly.express as px
from dash import dcc

from covigator.database.model import PrecomputedVariantsPerLineage, Variant


class LineageFigures(Figures):

    def get_lineages_plot(self, data_source: str, countries=None, lineages=None):
        logger.debug("Getting data on samples by country...")
        data = self.queries.get_accumulated_lineages_by_country(
            data_source=data_source, countries=countries, lineages=lineages)
        graph = dcc.Markdown("""**No data for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Prepare plot on samples by lineage...")
            lineages = list(data.sort_values("cumsum", ascending=False).lineage.unique())

            fig = make_subplots(rows=2, cols=1)

            fig1 = px.area(
                data, x="date", y="cumsum", color="lineage",
                category_orders={"lineage": lineages[::-1]},
                labels={"ratio_per_date": "% daily samples", "cumsum": "num. samples", "count": "increment"},
                hover_data=["ratio_per_date", "count"],
                color_discrete_sequence=px.colors.qualitative.Vivid)
            fig1.update_traces(line=dict(width=0.5), showlegend=False)

            fig2 = px.area(
                data, x="date", y="ratio_per_date", color="lineage",
                category_orders={"lineage": lineages[::-1]},
                labels={"ratio_per_date": "% daily samples", "cumsum": "num. samples", "count": "increment"},
                hover_data=["cumsum", "count"],
                color_discrete_sequence=px.colors.qualitative.Vivid)
            fig2.update_traces(line=dict(width=0.5))

            for trace in fig1["data"]:
                fig.append_trace(trace, row=1, col=1)
            for trace in fig2["data"]:
                fig.append_trace(trace, row=2, col=1)

            fig.update_layout(
                margin=MARGIN,
                template=TEMPLATE,
                legend={'traceorder': 'reversed', 'title': None},
                xaxis={'title': None},
                yaxis={
                    'title': 'num. samples'
                },
                yaxis2={
                    'tickformat': ',.0%',
                    'range': [0, 1],
                    'title': '% samples'
                },
                height=700,
            )
            fig.update_xaxes(showspikes=True)

            top_lineages = lineages[0: min(5, len(lineages))]
            top_lineages_and_cumsum = data[data.lineage.isin(top_lineages)][["lineage", "cumsum"]] \
                .groupby("lineage").max().reset_index().sort_values("cumsum", ascending=False)
            top_lineages_tooltip = list(
                top_lineages_and_cumsum.apply(lambda x: "{} ({})".format(x[0], int(x[1])), axis=1))

            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                **Samples by lineages**
                
                Top {} lineages: {}.
                """.format(len(top_lineages_tooltip),
                           ", ".join(top_lineages_tooltip)))
            ]
        return graph

    def discrete_background_color_bins(self, df, n_bins=5, columns='all', colors='Blues'):
        bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
        if columns == 'all':
            if 'id' in df:
                df_numeric_columns = df.select_dtypes('number').drop(['id'], axis=1)
            else:
                df_numeric_columns = df.select_dtypes('number')
        else:
            df_numeric_columns = df[columns]
        df_max = df_numeric_columns.max().max()
        df_min = df_numeric_columns.min().min()
        ranges = [
            ((df_max - df_min) * i) + df_min
            for i in bounds
        ]
        styles = []
        for i in range(1, len(bounds)):
            min_bound = ranges[i - 1]
            max_bound = ranges[i]
            backgroundColor = colorlover.scales[str(n_bins)]['seq'][colors][i - 1]
            color = 'white' if i > len(bounds) / 2. else 'inherit'

            for column in df_numeric_columns:
                styles.append({
                    'if': {
                        'filter_query': (
                                '{{{column}}} >= {min_bound}' +
                                (' && {{{column}}} < {max_bound}' if (i < len(bounds) - 1) else '')
                        ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                        'column_id': column
                    },
                    'backgroundColor': backgroundColor,
                    'color': color
                })

        return styles

    def get_lineages_variants_table(self, data_source: str, lineages: List[str] = None, countries: List[str] = None):

        result = None
        if lineages is not None and len(lineages) == 1:
            query = self.queries.session.query(
                PrecomputedVariantsPerLineage.variant_id,
                PrecomputedVariantsPerLineage.count_observations,
                #Variant.variant_type,
                Variant.gene_name,
                Variant.hgvs_p,
                Variant.annotation_highest_impact
            )\
                .filter(and_(
                    PrecomputedVariantsPerLineage.source == data_source,
                    PrecomputedVariantsPerLineage.lineage == lineages[0]))

            if countries is not None and len(countries) > 0:
                query = query.filter(PrecomputedVariantsPerLineage.country.in_(countries))

            query = query.join(Variant, PrecomputedVariantsPerLineage.variant_id == Variant.variant_id)\
                .order_by(Variant.position.asc())
            data = pd.read_sql(query.statement, self.queries.session.bind)

            # adds together the counts from multiple countries
            data = data.groupby(["variant_id", "gene_name", "hgvs_p", "annotation_highest_impact"]).sum().reset_index()

            styles_counts = self.discrete_background_color_bins(data, columns=["count_observations"])

            fig = dash_table.DataTable(
                id='lineages-variants-table',
                data=data.to_dict('records'),
                columns=[
                    {"name": ["Gene"], "id": "gene_name"},
                    {"name": ["DNA mutation"], "id": "variant_id"},
                    {"name": ["Protein change"], "id": "hgvs_p"},
                    {"name": ["Annotation"], "id": "annotation_highest_impact"},
                    {"name": ["Count observations"], "id": "count_observations"}
                ],
                style_data_conditional=STYLES_STRIPPED + styles_counts,
                style_as_list_view=True,
                style_header=STYLE_HEADER,
                style_cell=STYLE_CELL,
                css=[{'selector': '.dash-cell div.dash-cell-value',
                      'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'}],
                style_table={'overflowX': 'auto'},
                style_data={'whiteSpace': 'normal', 'height': 'auto'},
            )
            result = [fig, dcc.Markdown("**Mutations in lineage {lineage}**".format(lineage=lineages[0]))]

        return result
