from itertools import cycle
import colorlover
from dash import dash_table
import numpy as np
import re

import pandas as pd
from logzero import logger
from scipy.spatial.distance import squareform
from skbio.stats.distance import DissimilarityMatrix
from sklearn.cluster import OPTICS

from covigator import MISSENSE_VARIANT, DISRUPTIVE_INFRAME_DELETION, CONSERVATIVE_INFRAME_DELETION, \
    CONSERVATIVE_INFRAME_INSERTION, DISRUPTIVE_INFRAME_INSERTION
from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, TEMPLATE, MARGIN, STYLES_STRIPPED, STYLE_HEADER, \
    STYLE_CELL
import plotly
import plotly.express as px
import plotly.graph_objects as go
from dash import html
from dash import dcc
from covigator.database.model import Gene, Domain, DataSource

VARIANT_TOOLTIP = '<b>%{text}</b><br>' + 'Allele frequency: %{y:.5f}<br>' + 'Genomic Position: %{x}'
GENE_COLORS = list(reversed(plotly.express.colors.sequential.Tealgrn))
DOMAIN_COLORS = list(reversed(plotly.express.colors.sequential.Magenta))
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


class RecurrentMutationsFigures(Figures):

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

    def get_top_occurring_variants_plot(self, top, gene_name, domain, date_range_start, date_range_end, metric, source):
        data = self.queries.get_top_occurring_variants_precomputed(top, gene_name, domain, metric, source)

        fig = [
            dash_table.DataTable(id="top-occurring-variants-table"),
            dcc.Markdown("""**No mutations for the current selection**""")
        ]

        if data is not None and data.shape[0] > 0:
            # removes the columns from the months out of the range
            month_columns = [c for c in data.columns if MONTH_PATTERN.match(c)]
            included_month_colums = [c for c in month_columns if c >= date_range_start and c <= date_range_end]
            excluded_month_colums = [c for c in month_columns if c < date_range_start or c > date_range_end]
            data.drop(excluded_month_colums, axis=1, inplace=True)

            # set the styles of the cells
            styles_counts = self.discrete_background_color_bins(data, columns=included_month_colums)
            styles_total_count = self.discrete_background_color_bins(data, columns=["total"], colors="Reds")
            styles_frequency = self._get_table_style_by_af()
            month_columns = [{'name': ["", i], 'id': i} for i in data.columns if i.startswith("20")]
            month_columns[0]['name'][0] = 'Monthly counts' if metric == "count" else 'Monthly frequencies'

            fig = [
                dash_table.DataTable(
                    # this cannot be imported from the tab definition due to a circular import...
                    id="top-occurring-variants-table",
                    data=data.to_dict('records'),
                    sort_action='native',
                    columns=[
                                {"name": ["Variant", "Gene"], "id": "gene_name"},
                                {"name": ["", "DNA mutation"], "id": "dna_mutation"},
                                {"name": ["", "Protein mutation"], "id": "hgvs_p"},
                                {"name": ["", "Effect"], "id": "annotation"},
                                {"name": ["", "Frequency"], "id": "frequency"},
                                {"name": ["", "Count"], "id": "total"},
                            ] + month_columns,
                    style_data_conditional=STYLES_STRIPPED + styles_counts + styles_frequency + styles_total_count,
                    style_table={'overflowX': 'auto'},
                    style_as_list_view=True,
                    style_header=STYLE_HEADER,
                    style_cell=STYLE_CELL,
                    sort_by=[{"column_id": "frequency", "direction": "desc"}],
                    row_selectable='multi'
                ),
                html.Br(),
                html.Div(children=[
                    html.Button("Download CSV", id="btn_csv2"),
                    dcc.Download(id="download-dataframe-csv2"),
                    dcc.Store(id="memory2", data=data.to_dict('records'))]),
                html.Br(),
                dcc.Markdown("""
                    **Top occurring mutations table** 
                    *table shows the {} mutations{} with the highest frequency across all samples.*
                    *The counts and frequencies per month are only shown between {} and {}.*
                    *Selections in this table will be highlighted in the genome view and in the co-occurrence matrix.*
                    """.format(top, " in gene {}".format(gene_name) if gene_name else "",
                               date_range_start, date_range_end))
            ]

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

    def get_cooccurrence_heatmap(self, sparse_matrix, selected_variants, metric="jaccard", min_cooccurrences=5):
        data = self._get_variants_cooccurrence_matrix(data=sparse_matrix)
        graph = [dcc.Markdown("""**No co-occurrent mutations for the current selection**""")]
        if data is not None and data.shape[0] > 0:

            all_variants = data.variant_id_one.unique()
            values = np.array_split(data[metric], len(all_variants))
            texts = np.array_split(data.hgvs_tooltip, len(all_variants))
            if metric == "count":
                hovertemplate = '<b>%{text}</b><br>' + 'Counts: %{z}<br>' + 'Variant one: %{x}<br>' + 'Variant two: %{y}'
            elif metric == "frequency":
                hovertemplate = '<b>%{text}</b><br>' + 'Frequency: %{z:.3f}<br>' + 'Variant one: %{x}<br>' + 'Variant two: %{y}'
            elif metric == "jaccard":
                hovertemplate = '<b>%{text}</b><br>' + 'Jaccard index: %{z:.3f}<br>' + 'Variant one: %{x}<br>' + 'Variant two: %{y}'
            elif metric == "kappa":
                hovertemplate = '<b>%{text}</b><br>' + 'Kappa coefficient: %{z:.3f}<br>' + 'Variant one: %{x}<br>' + 'Variant two: %{y}'
            heatmap = go.Heatmap(
                z=values,
                x=all_variants,
                y=all_variants,
                colorscale="Oryel",
                hoverongaps=False,
                text=texts,
                hovertemplate=hovertemplate,
            )
            if selected_variants:
                selected_variant_ids = [v.get("dna_mutation") for v in selected_variants]
                selected_data = data
                selected_data[metric].where(
                    (selected_data.variant_id_one.isin(selected_variant_ids)) |
                    (selected_data.variant_id_two.isin(selected_variant_ids)), np.nan, inplace=True)
                values_selected = np.array_split(selected_data[metric], len(all_variants))
                texts_selected = np.array_split(selected_data.hgvs_tooltip, len(all_variants))

                heatmap_selected = go.Heatmap(
                    z=values_selected,
                    x=all_variants,
                    y=all_variants,
                    colorscale="Teal",
                    hoverongaps=False,
                    showscale=False,
                    text=texts_selected,
                    hovertemplate=hovertemplate,
                )
            layout = go.Layout(
                template=TEMPLATE,
                margin=MARGIN,
                height=700,
                yaxis=dict(visible=True, tickfont={"size": 10}, showgrid=False, showspikes=True, spikemode='toaxis',
                           spikethickness=2),
                xaxis=dict(tickangle=-45, tickfont={"size": 10}, showgrid=False, showspikes=True, spikemode='toaxis',
                           spikethickness=2),
            )
            traces = [heatmap]
            if selected_variants:
                traces.append(heatmap_selected)
            fig = go.Figure(data=traces, layout=layout)

            # the y index is reversed in plotly heatmap
            fig.update_yaxes(autorange="reversed")
            graph = [
                dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                dcc.Markdown("""
                                    ***Co-occurrence matrix*** *showing variant pairs co-occurring in at least {} samples (this value is configurable).*
                                    *The metric in the co-occurrence matrix can be chosen among counts, frequencies, Jaccard index or 
                                    Cohen's kappa coefficient. The Cohen's kappa coefficient introduces a correction to the Jaccard index for
                                    mutations with low occurrence.*
                                    *The diagonal contains the total counts or just 1.0 in the other metrics.*
                                    *The upper diagonal is not shown for clarity.*
                                    *Synonymous mutations are excluded.*
                                    *Different genomic mutations causing the same protein variant are not grouped.*
                                    """.format(min_cooccurrences))
            ]

        return html.Div(children=graph)

    def _get_variants_cooccurrence_matrix(self, data) -> pd.DataFrame:
        """
        Returns the full cooccurrence matrix of all non synonymous variants in a gene with at least
        min_occurrences occurrences.
        """
        # query for total samples required to calculate frequencies
        count_samples = self.queries.count_samples(source=DataSource.ENA.name)

        full_matrix = None
        if data.shape[0] > 0:
            # these are views of the original data
            annotations = data.loc[data.variant_id_one == data.variant_id_two,
                                   ["variant_id_one", "position", "reference", "alternate", "hgvs_p"]]
            tooltip = data.loc[:, ["variant_id_one", "variant_id_two", "hgvs_tooltip"]]
            diagonal = data.loc[
                data.variant_id_one == data.variant_id_two, ["variant_id_one", "variant_id_two", "count"]]
            sparse_matrix = data.loc[:, ["variant_id_one", "variant_id_two", "count"]]
            sparse_matrix["frequency"] = sparse_matrix["count"] / count_samples
            sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_one", how="left",
                                     suffixes=("", "_one"))
            sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_two", how="left",
                                     suffixes=("", "_two"))

            # calculate Jaccard index
            sparse_matrix["count_union"] = sparse_matrix["count_one"] + sparse_matrix["count_two"] - sparse_matrix[
                "count"]
            sparse_matrix["jaccard"] = sparse_matrix["count"] / sparse_matrix["count_union"]

            # calculate Cohen's kappa
            sparse_matrix["chance_agreement"] = np.exp(-sparse_matrix["count"])
            sparse_matrix["kappa"] = 1 - ((1 - sparse_matrix.jaccard) / (1 - sparse_matrix.chance_agreement))
            sparse_matrix["kappa"] = sparse_matrix["kappa"].transform(lambda k: k if k > 0 else 0)

            del sparse_matrix["count_union"]
            del sparse_matrix["count_one"]
            del sparse_matrix["count_two"]
            del sparse_matrix["chance_agreement"]

            # from the sparse matrix builds in memory the complete matrix
            all_variants = data.variant_id_one.unique()
            empty_full_matrix = pd.DataFrame(index=pd.MultiIndex.from_product(
                [all_variants, all_variants], names=["variant_id_one", "variant_id_two"])).reset_index()
            full_matrix = pd.merge(
                left=empty_full_matrix, right=sparse_matrix, on=["variant_id_one", "variant_id_two"], how='left')

            # add annotation on variant one
            full_matrix = pd.merge(left=full_matrix, right=annotations, on="variant_id_one", how='left') \
                .rename(columns={"position": "position_one",
                                 "reference": "reference_one",
                                 "alternate": "alternate_one",
                                 "hgvs_p": "hgvs_p_one"})

            # add annotations on variant two
            full_matrix = pd.merge(left=full_matrix, right=annotations,
                                   left_on="variant_id_two", right_on="variant_id_one", how='left') \
                .rename(columns={"variant_id_one_x": "variant_id_one",
                                 "position": "position_two",
                                 "reference": "reference_two",
                                 "alternate": "alternate_two",
                                 "hgvs_p": "hgvs_p_two"})

            # add tooltip
            full_matrix = pd.merge(left=full_matrix, right=tooltip, on=["variant_id_one", "variant_id_two"], how='left')
            # correct tooltip in diagonal
            full_matrix.loc[
                full_matrix.variant_id_one == full_matrix.variant_id_two, "hgvs_tooltip"] = full_matrix.hgvs_p_one

            # NOTE: transpose matrix manually as plotly transpose does not work with labels
            # the database return the upper diagonal, the lower is best for plots
            full_matrix.sort_values(["position_two", "reference_two", "alternate_two",
                                     "position_one", "reference_one", "alternate_one"], inplace=True)

            full_matrix = full_matrix.loc[:, ["variant_id_one", "variant_id_two", "count", "frequency", "jaccard",
                                              "kappa", "hgvs_tooltip"]]

        return full_matrix

    def get_variants_abundance_plot(self, source, bin_size=50):

        # reads genes and domains across the whole genome
        genes = self.queries.get_genes()
        domains = self.queries.get_domains()
        genes_and_domains = [(list(filter(lambda g: g.name == d.gene_name, genes))[0], d) for d in domains]

        # reads variants abundance
        variant_abundance = self.queries.get_variant_abundance_histogram(bin_size=bin_size, source=source)

        # reads conservation and bins it
        conservation = self.queries.get_conservation_table(bin_size=bin_size)

        # joins variant abundance and conservation
        data = variant_abundance.set_index("position_bin").join(conservation.set_index("position_bin"))
        data.reset_index(inplace=True)
        data.fillna(0, inplace=True)

        layout = go.Layout(
            template=TEMPLATE,
            margin=MARGIN,
            xaxis=dict(domain=[0, 1.0], tickformat=',d', hoverformat=',d', ticksuffix=" bp", ticks="outside",
                       visible=True, anchor="y7", showspikes=True, spikemode='across', spikethickness=2),
            yaxis=dict(domain=[0.9, 1.0], anchor='x7'),
            yaxis2=dict(domain=[0.6, 0.9], anchor='x7'),
            yaxis3=dict(domain=[0.45, 0.6], anchor='x7', visible=False),
            yaxis4=dict(domain=[0.3, 0.45], anchor='x7', visible=False),
            yaxis5=dict(domain=[0.15, 0.3], anchor='x7', visible=False),
            yaxis6=dict(domain=[0.05, 0.1], anchor='x7', visible=False),
            yaxis7=dict(domain=[0.0, 0.05], anchor='x7', visible=False),
            legend={'traceorder': 'normal'}
        )
        gene_colors = cycle(GENE_COLORS)
        domain_colors = cycle(DOMAIN_COLORS)
        gene_traces = [
            self._get_gene_trace(g, start=g.start, end=g.end, color=c, yaxis='y6') for g, c in zip(genes, gene_colors)]
        domain_traces = [self._get_domain_trace(color=c, domain=d, gene=g, yaxis='y7')
                         for (g, d), c in zip(genes_and_domains, domain_colors)]
        conservation_traces = self._get_conservation_traces(
            conservation, xaxis='x', yaxis1='y3', yaxis2='y4', yaxis3='y5')
        mean_unique_variants_per_bin = data.count_unique_variants.mean()
        variant_counts_traces = [
            go.Scatter(x=data.position_bin, y=data.count_variant_observations,
                       name="All variants", text="All variants", showlegend=False,
                       line_color=plotly.express.colors.sequential.Blues[-2], line_width=1),
            go.Scatter(x=data.position_bin,
                       y=[mean_unique_variants_per_bin for _ in range(data.shape[0])],
                       yaxis='y2', name="Mean unique variants", text="Mean unique variants",
                       line_width=1, showlegend=False, line_color=plotly.express.colors.sequential.Blues[-3]),
            go.Scatter(x=data.position_bin, y=data.count_unique_variants, yaxis='y2',
                       name="Unique variants", text="Unique variants", showlegend=False, fill='tonexty',
                       line_color=plotly.express.colors.sequential.Blues[-4], line_width=1)
        ]

        fig = go.Figure(data=variant_counts_traces + conservation_traces + gene_traces + domain_traces, layout=layout)

        # add track names
        fig.add_annotation(x=0.98, y=1.1, xref="x domain", yref="y domain", text="All variants",
                           showarrow=False, yshift=10)
        fig.add_annotation(x=0.98, y=0.9, xref="x domain", yref="y2 domain", text="Unique variants",
                           showarrow=False, yshift=10)
        fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y3 domain", text="Conservation SARS-CoV-2",
                           showarrow=False, yshift=10)
        fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y4 domain", text="Conservation SARS-like betaCoV",
                           showarrow=False, yshift=10)
        fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y5 domain", text="Conservation vertebrate CoV",
                           showarrow=False, yshift=10)

        return [
            dcc.Graph(
                figure=fig,
                config=PLOTLY_CONFIG
            ),
            dcc.Markdown("""
                ***Genome view*** *representing the abundance of mutations and ConsHMM (Arneson, 2019) conservation 
                using a bin size of {} bp. Synonymous mutations are included.*

                *The first track shows the count of mutations across the genome including repetitions. *
                *The second track shows the count of unique mutations across the genome, the horizontal line represents 
                the average number of unique mutations per bin and thus distinguishes regions with over and under the 
                average number of unique mutations.*
                *The third, fourth and fifth tracks represent the conservation as reported by ConsHMM within 
                SARS-CoV-2, among SARS-like betaCoV and among vertebrate CoV. Correlation between distribution of 
                unique mutations and conservation within Sars-CoV-2, among SARS-like betacoronavirus and among 
                vertebrates CoV respectively: {}, {}, {}.*
                *Genes and Pfam domains are represented in tones of red and purple respectively.*

                *Conservation data source: https://github.com/ernstlab/ConsHMM_CoV*

                *Arneson A, Ernst J. Systematic discovery of conservation states for single-nucleotide annotation of the 
                human genome. Communications Biology, 248, 2019. doi: https://doi.org/10.1038/s42003-019-0488-1*
                """.format(bin_size,
                           round(np.corrcoef(data.conservation, data.count_unique_variants)[0][1], 5),
                           round(np.corrcoef(data.conservation_sarbecovirus, data.count_unique_variants)[0][1], 5),
                           round(np.corrcoef(data.conservation_vertebrates, data.count_unique_variants)[0][1], 5)
                           ))]

    def get_variants_plot(self, gene_name, domain_name, selected_variants, bin_size, source):

        # reads gene annotations
        logger.debug("Getting genes and domains...")
        assert gene_name is not None or domain_name is not None, "Either gene or domain need to be provided"
        
        if domain_name is None:
            gene = self.queries.get_gene(gene_name)
            domains = self.queries.get_domains_by_gene(gene_name)
            start = gene.start
            end = gene.end
        else:
            domain = self.queries.get_domain(domain_name=domain_name)
            gene = self.queries.get_gene(domain.gene_name)
            domains = [domain]
            start = gene.start + (domain.start * 3)
            end = gene.start + (domain.end * 3)

        # reads variants
        logger.debug("Getting mutations...")
        variants = self.queries.get_non_synonymous_variants_by_region(start=start, end=end, source=source)

        # reads conservation and bins it
        logger.debug("Getting conservation...")
        conservation = self.queries.get_conservation_table(start=start, end=end, bin_size=bin_size)

        if variants.shape[0] > 0:
            # reads total number of samples and calculates frequencies
            logger.debug("Getting sample count...")
            count_samples = self.queries.count_samples(source=source)
            variants["af"] = variants.count_occurrences / count_samples
            variants["log_af"] = variants.af.transform(lambda x: np.log(x + 1))
            variants["log_count"] = variants.count_occurrences.transform(lambda x: np.log(x))
            variants.annotation_highest_impact = variants.annotation_highest_impact.transform(lambda a: a.split("&")[0])

            main_xaxis = 'x'

            variants_traces = []
            missense_variants = variants[variants.annotation_highest_impact == MISSENSE_VARIANT]
            if missense_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    missense_variants, name="missense variants", symbol=MISSENSE_VARIANT_SYMBOL, xaxis=main_xaxis))

            deletion_variants = variants[variants.annotation_highest_impact.isin(
                [DISRUPTIVE_INFRAME_DELETION, CONSERVATIVE_INFRAME_DELETION])]
            if deletion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    deletion_variants, name="inframe deletions", symbol=DELETION_SYMBOL, xaxis=main_xaxis))

            insertion_variants = variants[variants.annotation_highest_impact.isin(
                [DISRUPTIVE_INFRAME_INSERTION, CONSERVATIVE_INFRAME_INSERTION])]
            if insertion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    insertion_variants, name="inframe insertions", symbol=INSERTION_SYMBOL, xaxis=main_xaxis))

            other_variants = variants[~variants.annotation_highest_impact.isin([
                MISSENSE_VARIANT, DISRUPTIVE_INFRAME_DELETION, CONSERVATIVE_INFRAME_DELETION,
                DISRUPTIVE_INFRAME_INSERTION, CONSERVATIVE_INFRAME_DELETION])]
            if other_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    other_variants, name="other variants", symbol=OTHER_VARIANT_SYMBOL, xaxis=main_xaxis))

            selected_variants_trace = None
            if selected_variants:
                selected_variants_trace = go.Scatter(
                    x=[int(v.get("dna_mutation").split(":")[0]) for v in selected_variants],
                    y=[v.get("frequency") for v in selected_variants],
                    name="selected variants",
                    mode='markers',
                    marker=dict(
                        symbol="circle",
                        color="blue",
                        size=10,
                        showscale=False
                    ),
                    xaxis=main_xaxis,
                    showlegend=True,
                    text=["{} ({})".format(v.get("hgvs_p"), v.get("annotation")) for v in selected_variants],
                    hovertemplate=VARIANT_TOOLTIP
                )

            domain_colors = cycle(DOMAIN_COLORS)
            gene_trace = self._get_gene_trace(
                gene, start=start, end=end, color=plotly.express.colors.sequential.Tealgrn[-1], yaxis='y5', xaxis=main_xaxis)
            domain_traces = [self._get_domain_trace(
                color=c, gene=gene, domain=d, yaxis='y6', xaxis=main_xaxis, showlegend=True)
                for d, c in zip(domains, domain_colors)]
            conservation_traces = self._get_conservation_traces(
                conservation, main_xaxis, yaxis1='y2', yaxis2='y3', yaxis3='y4')

            data = variants_traces + conservation_traces + [gene_trace] + domain_traces
            if selected_variants_trace is not None:
                data.append(selected_variants_trace)

            layout = go.Layout(
                template=TEMPLATE,
                margin=go.layout.Margin(l=0, r=0, b=0, t=20),
                xaxis=dict(tickformat=',d', hoverformat=',d', ticksuffix=" bp", ticks="outside",
                           showspikes=True, spikemode='across', spikethickness=1, anchor='y6'),
                yaxis=dict(title='Allele frequency', type='log', domain=[0.4, 1.0], anchor=main_xaxis),
                yaxis2=dict(domain=[0.3, 0.4], visible=False, anchor=main_xaxis),
                yaxis3=dict(domain=[0.2, 0.3], visible=False, anchor=main_xaxis),
                yaxis4=dict(domain=[0.1, 0.2], visible=False, anchor=main_xaxis),
                yaxis5=dict(domain=[0.05, 0.1], visible=False, anchor=main_xaxis),
                yaxis6=dict(domain=[0.0, 0.05], visible=False, anchor=main_xaxis)
            )

            fig = go.Figure(data=data, layout=layout)

            fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y2 domain", text="Conservation SARS-CoV-2",
                               showarrow=False, yshift=10)
            fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y3 domain", text="Conservation SARS-like betaCoV",
                               showarrow=False, yshift=10)
            fig.add_annotation(x=0.98, y=1.0, xref="x domain", yref="y4 domain", text="Conservation vertebrate CoV",
                               showarrow=False, yshift=10)

            # add horizontal lines on the frequency boundaries
            fig.add_hline(y=0.1, line_width=1, line_dash="dash", line_color=VERY_COMMON_VARIANTS_COLOR)
            fig.add_hline(y=0.01, line_width=1, line_dash="dash", line_color=COMMON_VARIANTS_COLOR)
            fig.add_hline(y=0.001, line_width=1, line_dash="dash", line_color=RARE_VARIANTS_COLOR)

            return [dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
                    dcc.Markdown("""
                        ***Gene view*** *representing each variant with its frequency in the population and 
                        ConsHMM (Arneson, 2019) conservation using a bin size of {bin_size} bp. Synonymous mutations 
                        and mutations occurring in a single sample are excluded.*

                        *The scatter plot shows non synonymous mutations occurring in at least two samples on {region}. 
                        The x-axis shows the genomic coordinates and the y-axis shows the allele frequency.*
                        *The category of "other variants" includes frameshift indels, stop codon gain and lost and 
                        start lost mutations.*
                        *The mutations are colored according to their frequency as rare mutations (< 0.1 %), 
                        low frequency mutations (>= 0.1% and < 1%), common mutations (>= 1% and < 10%) and very common 
                        mutations (>= 10%). The second, third and fourth tracks in grey represent the conservation as 
                        reported by ConsHMM within SARS-CoV-2, among SARS-like betaCoV and among vertebrate CoV.*
                        *Genes and Pfam domains are represented in tones of red and purple respectively.*

                        *Conservation data source: https://github.com/ernstlab/ConsHMM_CoV*

                        *Arneson A, Ernst J. Systematic discovery of conservation states for single-nucleotide annotation of the 
                        human genome. Communications Biology, 248, 2019. doi: https://doi.org/10.1038/s42003-019-0488-1*
                        """.format(
                        bin_size=bin_size,
                        region="gene {}".format(gene_name) if domain_name is not None else
                                "domain {}: {}".format(gene_name, domain_name)))]
        else:
            return dcc.Markdown("""**No mutations for the current selection**""")

    def _get_conservation_traces(self, conservation, xaxis, yaxis1, yaxis2, yaxis3):
        return [
            go.Scatter(x=conservation.position_bin, y=conservation.conservation, yaxis=yaxis1, xaxis=xaxis,
                       text="Conservation SARS-CoV-2", textposition="top right", showlegend=False,
                       fill='tozeroy', line_color="grey", line_width=1),
            go.Scatter(x=conservation.position_bin, y=conservation.conservation_sarbecovirus, yaxis=yaxis2,
                       xaxis=xaxis, text="Conservation SARS-like betacoronavirus", textposition="top right",
                       showlegend=False, fill='tozeroy', line_color="grey", line_width=1),
            go.Scatter(x=conservation.position_bin, y=conservation.conservation_vertebrates, yaxis=yaxis3,
                       xaxis=xaxis, text="Conservation vertebrates", textposition="top right",
                       showlegend=False, fill='tozeroy', line_color="grey", line_width=1)
        ]

    @staticmethod
    def _get_domain_trace(color: str, domain: Domain, gene: Gene, yaxis='y', xaxis='x', showlegend=False):
        domain_start = gene.start + (domain.start * 3)  # domain coordinates are in the protein space
        domain_end = gene.start + (domain.end * 3)      # domain coordinates are in the protein space
        domain_name = domain.name
        return go.Scatter(
            mode='lines',
            x=[domain_start, domain_end, domain_end, domain_start],
            y=[0, 0, 1, 1],
            name=domain_name,
            fill="toself",
            fillcolor=color,
            text="<b>{} {}</b>: {}-{}".format(gene.name, domain_name, domain_start, domain_end),
            hoverinfo="text",
            line=dict(width=0),
            yaxis=yaxis,
            xaxis=xaxis,
            legendgroup='domains',
            showlegend=showlegend
        )

    @staticmethod
    def _get_gene_trace(gene: Gene, start, end, color, yaxis="y", xaxis='x'):
        # start and end coordinates in some case will come from the domain
        return go.Scatter(
            mode='lines',
            x=[start, end, end, start],
            y=[0, 0, 1, 1],
            name=gene.name,
            fill="toself",
            fillcolor=color,
            text="<b>{}</b>: {}-{}".format(gene.name, gene.start, gene.end),
            hoverinfo="text",
            line=dict(width=0),
            yaxis=yaxis,
            xaxis=xaxis,
            legendgroup='gene'
        )

    def _get_variants_scatter(self, variants, name, symbol, xaxis='x'):
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
            xaxis=xaxis,
            showlegend=True,
            text=variants[["hgvs_p", "annotation_highest_impact"]].apply(lambda x: "{} ({})".format(x[0], x[1]), axis=1),
            hovertemplate=VARIANT_TOOLTIP
        )

    def get_variants_clustering(self, sparse_matrix, min_cooccurrence, min_samples):

        data = self._run_clustering(sparse_matrix=sparse_matrix, min_samples=min_samples)

        tables = []
        if data is not None:
            for c in data.cluster.unique():
                tables.append(dash_table.DataTable(
                    id="cluster{}-variants-table".format(c),
                    data=data[data.cluster == c].to_dict('records'),
                    columns=[
                        {"name": [
                            "Cluster {} (mean Jaccard={})".format(c, data[data.cluster == c].cluster_jaccard_mean.iloc[0]),
                            "DNA mutation"], "id": "variant_id"},
                        {"name": ["", "Protein mutation"], "id": "hgvs_p"},
                    ],
                    fixed_columns={'headers': True, 'data': 1},
                    style_table={'overflowX': 'auto'},
                    style_cell={'minWidth': '50px', 'width': '50px', 'maxWidth': '50px'},
                    style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'fontWeight': 'bold'
                    },
                    sort_by=[{"column_id": "variant_id", "direction": "asc"}],
                ))
                tables.append(html.Br())
            tables.append(html.Div(children=[
                html.Button("Download CSV", id="btn_csv"),
                dcc.Download(id="download-dataframe-csv"),
                dcc.Store(id="memory", data=data.to_dict('records'))]))
            tables.append(html.Br()),
            tables.append(dcc.Markdown("""
            ***Co-occurrence clustering*** *shows the resulting clusters from the
            co-occurrence matrix with the Jaccard index corrected with the Cohen's kappa coefficient. 
            The co-occurrence matrix is built taking into account only mutations with at least {} pairwise 
            co-occurrences and if a gene is provided only mutations within that gene.  
            Clustering is performed on the co-occurrence matrix using OPTICS. 
            The mimimum number of neighbours to call
            a cluster is {}.
            Variants selected in the top occurrent mutations table are highlighted with a greater size in the plot.
            *

            *Ankerst et al. “OPTICS: ordering points to identify the clustering structure.” ACM SIGMOD Record 28, no. 2 (1999): 49-60.*
            """.format(min_cooccurrence, min_samples)))

        return html.Div(children=tables)

    def _run_clustering(self, sparse_matrix, min_samples) -> pd.DataFrame:

        data = None

        if sparse_matrix is not None and sparse_matrix.shape[0] > 0:
            diagonal = sparse_matrix.loc[sparse_matrix.variant_id_one == sparse_matrix.variant_id_two,
                                         ["variant_id_one", "variant_id_two", "count"]]
            sparse_matrix_with_diagonal = pd.merge(
                left=sparse_matrix, right=diagonal, on="variant_id_one", how="left", suffixes=("", "_one"))
            sparse_matrix_with_diagonal = pd.merge(
                left=sparse_matrix_with_diagonal, right=diagonal, on="variant_id_two", how="left", suffixes=("", "_two"))

            # calculate Jaccard index
            sparse_matrix_with_diagonal["count_union"] = sparse_matrix_with_diagonal["count_one"] + \
                                                         sparse_matrix_with_diagonal["count_two"] - \
                                                         sparse_matrix_with_diagonal["count"]
            sparse_matrix_with_diagonal["jaccard_similarity"] = sparse_matrix_with_diagonal["count"] / \
                                                                sparse_matrix_with_diagonal.count_union
            sparse_matrix_with_diagonal["jaccard_dissimilarity"] = 1 - sparse_matrix_with_diagonal.jaccard_similarity

            # calculate Cohen's kappa
            sparse_matrix_with_diagonal["chance_agreement"] = np.exp(-sparse_matrix_with_diagonal["count"])
            sparse_matrix_with_diagonal["kappa"] = 1 - ((1 - sparse_matrix_with_diagonal.jaccard_similarity) / (
                    1 - sparse_matrix_with_diagonal.chance_agreement))
            sparse_matrix_with_diagonal["kappa"] = sparse_matrix_with_diagonal["kappa"].transform(
                lambda k: k if k > 0 else 0)
            sparse_matrix_with_diagonal["kappa_dissimilarity"] = 1 - sparse_matrix_with_diagonal.kappa

            dissimilarity_metric = "kappa_dissimilarity"

            # build upper diagonal matrix
            all_variants = sparse_matrix_with_diagonal.variant_id_one.unique()
            empty_full_matrix = pd.DataFrame(index=pd.MultiIndex.from_product(
                [all_variants, all_variants], names=["variant_id_one", "variant_id_two"])).reset_index()
            upper_diagonal_matrix = pd.merge(
                # gets only the inferior matrix without the diagnonal
                left=empty_full_matrix.loc[empty_full_matrix.variant_id_one < empty_full_matrix.variant_id_two, :],
                right=sparse_matrix_with_diagonal.loc[:, ["variant_id_one", "variant_id_two", dissimilarity_metric]],
                on=["variant_id_one", "variant_id_two"], how='left')
            upper_diagonal_matrix.fillna(1.0, inplace=True)
            upper_diagonal_matrix.sort_values(by=["variant_id_one", "variant_id_two"], inplace=True)

            logger.debug("Building square distance matrix...")
            distance_matrix = squareform(upper_diagonal_matrix[dissimilarity_metric])
            # this ensures the order of variants ids is coherent with the non redundant form of the distance matrix
            ids = np.array([list(upper_diagonal_matrix.variant_id_one[0])[0]] + \
                           list(upper_diagonal_matrix.variant_id_two[0:len(upper_diagonal_matrix.variant_id_two.unique())]))
            distance_matrix_with_ids = DissimilarityMatrix(data=distance_matrix, ids=ids)

            logger.debug("Clustering...")
            clusters = OPTICS(min_samples=min_samples, max_eps=1.4).fit_predict(distance_matrix_with_ids.data)

            logger.debug("Building clustering dataframe...")
            data = pd.DataFrame({'variant_id': distance_matrix_with_ids.ids, 'cluster': clusters})

            logger.debug("Annotate with HGVS.p ...")
            annotations = pd.concat([
                sparse_matrix.loc[:, ["variant_id_one", "hgvs_p_one"]].rename(
                    columns={"variant_id_one": "variant_id", "hgvs_p_one": "hgvs_p"}),
                sparse_matrix.loc[:, ["variant_id_two", "hgvs_p_two"]].rename(
                    columns={"variant_id_two": "variant_id", "hgvs_p_two": "hgvs_p"})])
            data = pd.merge(left=data, right=annotations, on="variant_id", how="left")

            logger.debug("Annotate with cluster mean Jaccard index...")
            data["cluster_jaccard_mean"] = 1.0
            for c in data.cluster.unique():
                variants_in_cluster = data[data.cluster == c].variant_id.unique()
                data.cluster_jaccard_mean = np.where(data.cluster == c, sparse_matrix_with_diagonal[
                    (sparse_matrix.variant_id_one.isin(variants_in_cluster)) &
                    (sparse_matrix.variant_id_two.isin(variants_in_cluster)) &
                    (sparse_matrix.variant_id_one != sparse_matrix.variant_id_two)
                    ].jaccard_dissimilarity.mean(), data.cluster_jaccard_mean)

            data.cluster_jaccard_mean = data.cluster_jaccard_mean.transform(lambda x: round(x, 3))

            data = data.set_index("variant_id").drop_duplicates().reset_index().sort_values("cluster")
            data = data[data.cluster != -1].loc[:, ["cluster", "variant_id", "hgvs_p", "cluster_jaccard_mean"]]

        return data
