import random
import colorlover
import dash_table
import numpy as np
import plotly
import plotly.express as px
import plotly.graph_objects as go
import dash_core_components as dcc
import json

from six.moves.urllib import request as urlreq
import dash_bio as dashbio
import re
from covigator.database.queries import Queries

OTHER_VARIANT_SYMBOL = "x"

INSERTION_SYMBOL = "triangle-up"

DELETION_SYMBOL = "triangle-down"

MISSENSE_VARIANT_SYMBOL = "circle"

VERY_COMMON_VARIANTS_COLOR = plotly.express.colors.sequential.Reds[-1]
COMMON_VARIANTS_COLOR = plotly.express.colors.sequential.Reds[-3]
RARE_VARIANTS_COLOR = plotly.express.colors.sequential.Reds[-7]
COMMON_VARIANTS_THRESHOLD = 0.1
LOW_FREQUENCY_VARIANTS_THRESHOLD = 0.01
LOW_FREQUENCY_VARIANTS_COLOR = plotly.express.colors.sequential.Reds[-5]
RARE_VARIANTS_THRESHOLD = 0.001
MONTH_PATTERN = re.compile('[0-9]{4}-[0-9]{2}')


class Figures:

    def __init__(self, queries: Queries):
        self.queries = queries

    def get_accumulated_samples_by_country_plot(self):

        data = self.queries.get_accumulated_samples_by_country()
        fig = None
        if data is not None:
            fig = px.area(data, x="date", y="cumsum", color="country",
                          category_orders={
                              "country": list(data.sort_values("cumsum", ascending=False).country.unique())[::-1]},
                          labels={"cumsum": "num. samples", "count": "increment"},
                          title="Accumulated samples per country",
                          hover_data=["count"],
                          color_discrete_sequence=random.shuffle(px.colors.qualitative.Dark24))
            fig.update_layout(
                legend={'traceorder': 'reversed'},
                xaxis={'title': None},
                yaxis={'dtick': 2000}
            )
        return fig

    def _get_color_by_af(self, af):
        color = None
        if af < RARE_VARIANTS_THRESHOLD:
            color = RARE_VARIANTS_COLOR
        elif RARE_VARIANTS_THRESHOLD <= af < LOW_FREQUENCY_VARIANTS_THRESHOLD:
            color = LOW_FREQUENCY_VARIANTS_COLOR
        elif LOW_FREQUENCY_VARIANTS_THRESHOLD <= af < COMMON_VARIANTS_THRESHOLD:
            color = COMMON_VARIANTS_COLOR
        elif af >= COMMON_VARIANTS_THRESHOLD:
            color = VERY_COMMON_VARIANTS_COLOR
        return color

    def _get_table_style_by_af(self):
        return [
            {
                'if': {
                    'filter_query': '{{frequency}} >= 0 && {{frequency}} < {}'.format(RARE_VARIANTS_THRESHOLD),
                    'column_id': "frequency"
                },
                'backgroundColor': RARE_VARIANTS_COLOR,
                'color': 'inherit'
            },
            {
                'if': {
                    'filter_query': '{{frequency}} >= {} && {{frequency}} < {}'.format(
                        RARE_VARIANTS_THRESHOLD, LOW_FREQUENCY_VARIANTS_THRESHOLD),
                    'column_id': "frequency"
                },
                'backgroundColor': LOW_FREQUENCY_VARIANTS_COLOR,
                'color': 'inherit'
            },
            {
                'if': {
                    'filter_query': '{{frequency}} >= {} && {{frequency}} < {}'.format(
                        LOW_FREQUENCY_VARIANTS_THRESHOLD, COMMON_VARIANTS_THRESHOLD),
                    'column_id': "frequency"
                },
                'backgroundColor': COMMON_VARIANTS_COLOR,
                'color': 'white'
            },
            {
                'if': {
                    'filter_query': '{{frequency}} >= {}'.format(COMMON_VARIANTS_THRESHOLD),
                    'column_id': "frequency"
                },
                'backgroundColor': VERY_COMMON_VARIANTS_COLOR,
                'color': 'white'
            }
        ]

    def _get_variants_scatter(self, variants, name, symbol):
        return go.Scatter(
            x=variants.position,
            y=variants.af,
            name=name,
            mode='markers',
            # opacity=0.5,
            marker=dict(
                symbol=symbol,
                color=variants.af.transform(self._get_color_by_af),
                showscale=False
            ),
            xaxis='x1',
            showlegend=True,
            text=variants[["hgvs_p", "annotation"]].apply(lambda x: "{} ({})".format(x[0], x[1]), axis=1),
            hovertemplate='<b>%{text}</b><br>' +
                          'Allele frequency: %{y:.2f}<br>' +
                          'Genomic Position: %{x}'
        )

    def get_variants_plot(self, gene_name):

        # reads gene annotations
        gene = self.queries.get_gene(gene_name)
        gene_start = int(gene.data.get("start"))
        gene_end = int(gene.data.get("end"))
        protein_features = gene.data.get("transcripts", [])[0].get("translations", [])[0].get("protein_features")
        pfam_protein_features = [f for f in protein_features if f.get("dbname") == "Pfam"]

        # reads variants
        variants = self.queries.get_non_synonymous_variants_by_gene(gene_name)

        fig = dcc.Markdown("""**No variants for the current selection**""")
        if variants.shape[0] > 0:
            # reads total number of samples and calculates frequencies
            count_samples = self.queries.count_ena_samples_loaded()
            variants["af"] = variants.count_occurrences / count_samples
            variants["log_af"] = variants.af.transform(lambda x: np.log(x + 1))
            variants["log_count"] = variants.count_occurrences.transform(lambda x: np.log(x))
            # TODO: do something in the data ingestion about multiple annotations on the same variant
            variants.annotation = variants.annotation.transform(lambda a: a.split("&")[0])

            variants_traces = []
            missense_variants = variants[variants.annotation == "missense_variant"]
            if missense_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    missense_variants, name="missense variants", symbol=MISSENSE_VARIANT_SYMBOL))

            deletion_variants = variants[variants.annotation.isin(
                ["disruptive_inframe_deletion", "conservative_inframe_deletion"])]
            if deletion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    deletion_variants, name="inframe deletions", symbol=DELETION_SYMBOL))

            insertion_variants = variants[variants.annotation.isin(
                ["disruptive_inframe_insertion", "conservative_inframe_insertion"])]
            if insertion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    insertion_variants, name="inframe insertions", symbol=INSERTION_SYMBOL))

            other_variants = variants[~variants.annotation.isin([
                "missense_variant", "disruptive_inframe_deletion", "conservative_inframe_deletion",
                "disruptive_inframe_insertion", "conservative_inframe_insertion"])]
            if other_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    other_variants, name="other variants", symbol=OTHER_VARIANT_SYMBOL))

            gene_trace = go.Scatter(
                mode='lines',
                x=[gene_start, gene_end, gene_end, gene_start],
                y=[0.4, 0.4, 0.6, 0.6],
                name=gene_name,
                fill="toself",
                fillcolor="grey",
                hovertext=gene_name,
                hoveron="fills",
                line=dict(width=0),
                yaxis='y2',
                xaxis='x1',
                legendgroup='gene'
            )
            domain_traces = []
            for d, c in zip(pfam_protein_features, plotly.express.colors.qualitative.Plotly[0:len(pfam_protein_features)]):
                domain_start = gene_start + int(d["start"])
                domain_end = gene_start + int(d["end"])
                domain_name = d.get('description')
                domain_traces.append(go.Scatter(
                    mode='lines',
                    x=[domain_start, domain_end, domain_end, domain_start],
                    y=[0, 0, 1, 1],
                    name=domain_name,
                    fill="toself",
                    fillcolor=c,
                    hovertext=domain_name,
                    hoveron="fills",
                    line=dict(width=0),
                    yaxis='y2',
                    xaxis='x1',
                    legendgroup='domains'
                ))

            data = variants_traces + [gene_trace] + domain_traces
            layout = go.Layout(
                template="plotly_white",
                xaxis=dict(
                    domain=[0, 1.0],
                    tickformat=',d',
                    hoverformat=',d',
                    visible=False
                ),
                xaxis2=dict(
                    title='Genomic position',
                    tickformat=',d',
                    hoverformat=',d',
                    domain=[0, 1.0],
                    anchor='y2'
                ),
                yaxis=dict(
                    title='Allele frequency',
                    type='log',
                    domain=[0.1, 1.0],
                    anchor='x2'
                ),
                yaxis2=dict(
                    domain=[0.0, 0.1],
                    visible=False,
                    anchor='x2'
                ),
                margin=go.layout.Margin(l=0, r=0, b=0, t=0)
            )
            fig = go.Figure(data=data, layout=layout)

        return fig

    def get_top_occurring_variants_plot(self, top, gene_name, date_range_start, date_range_end):
        data = self.queries.get_top_occurring_variants(top, gene_name)
        fig = dcc.Markdown("""**No variants for the current selection**""")
        if data is not None and data.shape[0] > 0:

            # removes the columns from the months out of the range
            month_columns = [c for c in data.columns if MONTH_PATTERN.match(c)]
            included_month_colums = [c for c in month_columns if c >= date_range_start and c <= date_range_end]
            excluded_month_colums = [c for c in month_columns if c < date_range_start or c > date_range_end]
            data.drop(excluded_month_colums, axis=1, inplace=True)

            # set the styles of the cells
            styles_counts = self.discrete_background_color_bins(data, columns=included_month_colums)
            styles_frequency = self._get_table_style_by_af()
            styles_striped = [{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }]

            month_columns = [{'name': ["", i], 'id': i} for i in data.columns if i.startswith("20")]
            month_columns[0]['name'][0] = 'Monthly count'

            fig = dash_table.DataTable(
                    data=data.to_dict('records'),
                    sort_action='native',
                    columns=[
                        {"name": ["Variant", "Gene"], "id": "gene_name"},
                        {"name": ["", "DNA mutation"], "id": "dna_mutation"},
                        {"name": ["", "Protein mutation"], "id": "hgvs_p"},
                        {"name": ["", "Effect"], "id": "annotation"},
                        {"name": ["", "Frequency"], "id": "frequency"},
                    ] + month_columns,
                    style_data_conditional=styles_striped + styles_counts + styles_frequency,
                    style_cell_conditional=[
                        {
                            'if': {'column_id': c},
                            'textAlign': 'left'
                        } for c in ['gene_name', 'dna_mutation', 'hgvs_p', 'annotation']
                    ],
                    style_table={'overflowX': 'auto'},
                    style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'fontWeight': 'bold'
                    },
                    sort_by=[{"column_id": "frequency", "direction": "desc"}]
                )

        return fig

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

    def get_cooccurrence_heatmap(self, gene_name):
        data = self.queries.get_variants_cooccurrence_by_gene(gene_name=gene_name)
        fig = dcc.Markdown("""**No co-occurrent variants for the current selection**""")
        if data is not None and data.shape[0] > 0:
            data.sort_values(["variant_one", "variant_two"], inplace=True)
            all_variants = data.variant_one.unique()
            values = np.array_split(data["count"], len(all_variants))
            heatmap = go.Heatmap(
                z=values,
                x=all_variants,
                y=all_variants,
                colorscale="Greens",
                hoverongaps=False)
            layout = go.Layout(
                template="plotly_white",
                height=700,
                yaxis=dict(
                    #scaleanchor='x',
                    visible=False),
                # xaxis=dict(
                #     domain=[0, 1.0],
                #     tickformat=',d',
                #     hoverformat=',d',
                #     visible=False
                # ),
                # xaxis2=dict(
                #     title='Genomic position',
                #     tickformat=',d',
                #     hoverformat=',d',
                #     domain=[0, 1.0],
                #     anchor='y2'
                # ),
                # yaxis=dict(
                #     title='Allele frequency',
                #     type='log',
                #     domain=[0.1, 1.0],
                #     anchor='x2'
                # ),
                # yaxis2=dict(
                #     domain=[0.0, 0.1],
                #     visible=False,
                #     anchor='x2'
                # ),
                margin=go.layout.Margin(l=0, r=0, b=0, t=0)
            )
            fig = go.Figure(data=[heatmap], layout=layout)
        return fig
