from itertools import cycle
import colorlover
import dash_table
import numpy as np
import re
from covigator.dashboard.figures.figures import Figures, PLOTLY_CONFIG, TEMPLATE, MARGIN, STYLES_STRIPPED, STYLE_HEADER
import plotly
import plotly.express as px
import plotly.graph_objects as go
import dash_html_components as html
import dash_core_components as dcc
from covigator.database.model import Gene

VARIANT_TOOLTIP = '<b>%{text}</b><br>' + 'Allele frequency: %{y:.5f}<br>' + 'Genomic Position: %{x}'
GENE_COLORS = cycle(plotly.express.colors.sequential.Reds)
DOMAIN_COLORS = cycle(reversed(plotly.express.colors.sequential.Purples))
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


class VariantsFigures(Figures):

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

    def get_top_occurring_variants_plot(self, top, gene_name, date_range_start, date_range_end, metric, source):
        data = self.queries.get_top_occurring_variants_precomputed(top, gene_name, metric, source)
        fig = dcc.Markdown("""**No variants for the current selection**""")
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

            fig = dash_table.DataTable(
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
                style_cell_conditional=[
                    {
                        'if': {'column_id': c},
                        'textAlign': 'left'
                    } for c in ['gene_name', 'dna_mutation', 'hgvs_p', 'annotation']
                ],
                style_table={'overflowX': 'auto'},
                style_as_list_view=True,
                style_header=STYLE_HEADER,
                sort_by=[{"column_id": "frequency", "direction": "desc"}],
                row_selectable='multi'
            )

        return [
            fig,
            dcc.Markdown("""
                            **Top occurring variants** 
                            *table shows the {} variants{} with the highest frequency across all samples.*
                            *The counts and frequencies per month are only shown between {} and {}.*
                            *Selections in this table will be highlighted in the genome view and in the co-occurrence matrix.*
                            """.format(top, " in gene {}".format(gene_name) if gene_name else "",
                                       date_range_start, date_range_end))
        ]

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

    def get_cooccurrence_heatmap(self, gene_name, selected_variants, metric="jaccard", min_occurrences=5):
        data = self.queries.get_variants_cooccurrence_by_gene(gene_name=gene_name, min_cooccurrence=min_occurrences)
        graph = dcc.Markdown("""**No co-occurrent variants for the current selection**""")
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
            graph = dcc.Graph(figure=fig, config=PLOTLY_CONFIG)

        return [
            graph,
            dcc.Markdown("""
                        ***Co-occurrence matrix*** *showing variant pairs co-occurring in at least {} samples (this value is configurable).*
                        *The metric in the co-occurrence matrix can be chosen among counts, frequencies, Jaccard index or 
                        Cohen's kappa coefficient. The Cohen's kappa coefficient introduces a correction to the Jaccard index for
                        variants with low occurrence.*
                        *The diagonal contains the total counts or just 1.0 in the other metrics.*
                        *The upper diagonal is not shown for clarity.*
                        *Synonymous variants are excluded.*
                        *Different genomic variants causing the same protein variant are not grouped.*
                        """.format(min_occurrences, metric, " on gene {}".format(gene_name) if gene_name else ""))
        ]

    def get_variants_abundance_plot(self, bin_size=50, source=None):

        # reads genes and domains across the whole genome
        genes = self.queries.get_genes()
        domains = []
        for g in genes:
            domains.extend([(g, d) for d in g.get_pfam_domains()])

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

        gene_traces = [self._get_gene_trace(g, color=c, yaxis='y6') for g, c in zip(genes, GENE_COLORS)]
        domain_traces = [self._get_domain_trace(color=c, domain=d, gene=g, yaxis='y7')
                         for (g, d), c in zip(domains, DOMAIN_COLORS)]
        conservation_traces = self._get_conservation_traces(
            conservation, xaxis='x', yaxis1='y3', yaxis2='y4', yaxis3='y5')
        variant_counts_traces = [
            go.Scatter(x=data.position_bin, y=data.count_variant_observations,
                       name="All variants", text="All variants", showlegend=False,
                       line_color=plotly.express.colors.sequential.Blues[-2], line_width=1),
            go.Scatter(x=data.position_bin,
                       y=[data.count_unique_variants.mean() for _ in range(data.shape[0])],
                       yaxis='y2', name="Mean unique variants", text="Mean unique variants",
                       line_width=1,
                       showlegend=False, line_color=plotly.express.colors.sequential.Blues[-3]),
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
                ***Genome view*** *representing the abundance of variants and ConsHMM (Arneson, 2019) conservation 
                using a bin size of {} bp. Synonymous variants are included.*

                *The first track shows the count of variants across the genome including repetitions. *
                *The second track shows the count of unique variants across the genome, the horizontal line represents 
                the average number of unique variants per bin and thus distinguishes regions with over and under the 
                average number of unique variants.*
                *The third, fourth and fifth tracks represent the conservation as reported by ConsHMM within 
                SARS-CoV-2, among SARS-like betaCoV and among vertebrate CoV. Correlation between distribution of 
                unique variants and conservation within Sars-CoV-2, among SARS-like betacoronavirus and among 
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

    def get_variants_plot(self, gene_name, selected_variants, bin_size, source=None):

        # reads gene annotations
        gene = self.queries.get_gene(gene_name)
        pfam_protein_features = gene.get_pfam_domains()

        # reads variants
        variants = self.queries.get_non_synonymous_variants_by_region(start=gene.start, end=gene.end, source=source)

        # reads conservation and bins it
        conservation = self.queries.get_conservation_table(start=gene.start, end=gene.end, bin_size=bin_size)

        if variants.shape[0] > 0:
            # reads total number of samples and calculates frequencies
            count_samples = self.queries.count_samples(source=source)
            variants["af"] = variants.count_occurrences / count_samples
            variants["log_af"] = variants.af.transform(lambda x: np.log(x + 1))
            variants["log_count"] = variants.count_occurrences.transform(lambda x: np.log(x))
            # TODO: do something in the data ingestion about multiple annotations on the same variant
            variants.annotation = variants.annotation.transform(lambda a: a.split("&")[0])

            main_xaxis = 'x'

            variants_traces = []
            missense_variants = variants[variants.annotation == "missense_variant"]
            if missense_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    missense_variants, name="missense variants", symbol=MISSENSE_VARIANT_SYMBOL, xaxis=main_xaxis))

            deletion_variants = variants[variants.annotation.isin(
                ["disruptive_inframe_deletion", "conservative_inframe_deletion"])]
            if deletion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    deletion_variants, name="inframe deletions", symbol=DELETION_SYMBOL, xaxis=main_xaxis))

            insertion_variants = variants[variants.annotation.isin(
                ["disruptive_inframe_insertion", "conservative_inframe_insertion"])]
            if insertion_variants.shape[0] > 0:
                variants_traces.append(self._get_variants_scatter(
                    insertion_variants, name="inframe insertions", symbol=INSERTION_SYMBOL, xaxis=main_xaxis))

            other_variants = variants[~variants.annotation.isin([
                "missense_variant", "disruptive_inframe_deletion", "conservative_inframe_deletion",
                "disruptive_inframe_insertion", "conservative_inframe_insertion"])]
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

            gene_trace = self._get_gene_trace(
                gene, color=plotly.express.colors.sequential.Reds[1], yaxis='y5', xaxis=main_xaxis)
            domain_traces = [self._get_domain_trace(
                color=c, gene=gene, domain=d, yaxis='y6', xaxis=main_xaxis, showlegend=True)
                for d, c in zip(pfam_protein_features, DOMAIN_COLORS)]
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
                        ConsHMM (Arneson, 2019) conservation using a bin size of {} bp. Synonymous variants and variants 
                        occurring in a single sample are excluded.*

                        *The scatter plot shows non synonymous variants occurring in at least two samples on gene {}. 
                        The x-axis shows the genomic coordinates and the y-axis shows the allele frequency.*
                        *The category of "other variants" includes frameshift indels, stop codon gain and lost and 
                        start lost variants.*
                        *The variants are colored according to their frequency as rare variants (< 0.1 %), 
                        low frequency variants (>= 0.1% and < 1%), common variants (>= 1% and < 10%) and very common 
                        variants (>= 10%). The second, third and fourth tracks in grey represent the conservation as 
                        reported by ConsHMM within SARS-CoV-2, among SARS-like betaCoV and among vertebrate CoV.*
                        *Genes and Pfam domains are represented in tones of red and purple respectively.*

                        *Conservation data source: https://github.com/ernstlab/ConsHMM_CoV*

                        *Arneson A, Ernst J. Systematic discovery of conservation states for single-nucleotide annotation of the 
                        human genome. Communications Biology, 248, 2019. doi: https://doi.org/10.1038/s42003-019-0488-1*
                        """.format(bin_size, gene_name))]
        else:
            return dcc.Markdown("""**No variants for the current selection**""")

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
    def _get_domain_trace(color: str, domain: dict, gene: Gene, yaxis='y', xaxis='x', showlegend=False):
        domain_start = gene.start + int(domain["start"])
        domain_end = gene.start + int(domain["end"])
        domain_name = domain.get('description')
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
    def _get_gene_trace(gene: Gene, color, yaxis="y", xaxis='x'):
        return go.Scatter(
            mode='lines',
            x=[gene.start, gene.end, gene.end, gene.start],
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
            text=variants[["hgvs_p", "annotation"]].apply(lambda x: "{} ({})".format(x[0], x[1]), axis=1),
            hovertemplate=VARIANT_TOOLTIP
        )

    def get_variants_clustering(self, gene_name, selected_variants, min_cooccurrence, min_samples):
        data = self.queries.get_mds(
            gene_name=gene_name, min_cooccurrence=min_cooccurrence, min_samples=min_samples)

        tables = []
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

        return [
            html.Div(children=tables),
            dcc.Markdown("""
            ***Co-occurrence clustering*** *shows the resulting clusters from the
            co-occurrence matrix with the Jaccard index corrected with the Cohen's kappa coefficient. 
            The co-occurrence matrix is built taking into account only variants with at least {} pairwise 
            co-occurrences and if a gene is provided only variants within that gene.  
            Clustering is performed on the co-occurrence matrix using OPTICS. 
            The mimimum number of neighbours to call
            a cluster is {}.
            Variants selected in the top occurrent variants table are highlighted with a greater size in the plot.
            *

            *Ankerst et al. “OPTICS: ordering points to identify the clustering structure.” ACM SIGMOD Record 28, no. 2 (1999): 49-60.*
            """.format(min_cooccurrence, min_samples))]
