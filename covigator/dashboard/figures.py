import random
import colorlover
import dash_table
import numpy as np
import plotly.express as px
import dash_bio
import dash_core_components as dcc
import json
from six.moves.urllib import request as urlreq
import dash_bio as dashbio
import re
from covigator.database.queries import Queries


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

    def get_variants_plot(self, gene_name="S"):

        # reads gene annotations
        gene = self.queries.get_gene(gene_name)
        start = int(gene.data.get("start"))
        protein_features = gene.data.get("transcripts", [])[0].get("translations", [])[0].get("protein_features")
        domains = [{"name": f.get('description'), "coord": "{}-{}".format(
            start + int(f["start"]), start + int(f["end"]))} for f in protein_features if f.get("dbname") == "Pfam"]

        # reads variants
        variants = self.queries.get_non_synonymous_variants_by_gene(gene_name)

        fig = None
        if variants.shape[0] > 0:
            # reads total number of samples and calculates frequencies
            count_samples = self.queries.count_ena_samples_loaded()
            variants["af"] = variants.count_1 / count_samples
            variants["log_af"] = variants.af.transform(lambda x: np.log(x + 1))
            variants["log_count"] = variants.count_1.transform(lambda x: np.log(x))

            # TODO: do something in the data ingestion about multiple annotations on the same variant
            variants.annotation = variants.annotation.transform(lambda a: a.split("&")[0])

            mdata = {
                "x": list(variants.position.transform(lambda x: str(x))),
                "y": list(variants.log_count.transform(lambda x: round(x, 3))),
                "mutationGroups": list(variants.annotation),
                "domains": domains
            }

            fig = dash_bio.NeedlePlot(
                id='my-dashbio-needleplot',
                mutationData=mdata,
                rangeSlider=True,
                ylabel="log(count observations)",
                domainStyle={
                    'displayMinorDomains': False
                }
            )
        return fig

    def get_circos_plot(self):
        data = urlreq.urlopen(
            'https://raw.githubusercontent.com/plotly/dash-bio-docs-files/master/circos_graph_data.json').read()
        circos_graph_data = json.loads(data)

        layout = [{"id":"MN908947.3", "label": "MN908947.3", "color":"#996600", "len":29903}]

        return dashbio.Circos(
            layout=layout,
            tracks=[{
                'type': 'CHORDS',
                'data': circos_graph_data['chords']
            }],
            config={
                #'innerRadius': 40,
                #'outerRadius': 200
            }
        )

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
            styles_frequency = self.discrete_background_color_bins(data, columns=["frequency"], colors="Reds")
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
