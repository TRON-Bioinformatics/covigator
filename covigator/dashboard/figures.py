import os
import random
from itertools import cycle

import colorlover
import dash_table
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
import dash_core_components as dcc
import re

from covigator import ENV_COVIGATOR_STORAGE_FOLDER
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
            xaxis='x',
            showlegend=True,
            text=variants[["hgvs_p", "annotation"]].apply(lambda x: "{} ({})".format(x[0], x[1]), axis=1),
            hovertemplate='<b>%{text}</b><br>' +
                          'Allele frequency: %{y:.5f}<br>' +
                          'Genomic Position: %{x}'
        )

    def get_variants_plot(self, gene_name, selected_variants):

        # reads gene annotations
        gene = self.queries.get_gene(gene_name)
        gene_start = int(gene.data.get("start"))
        gene_end = int(gene.data.get("end"))
        pfam_protein_features = self.queries.get_pfam_domains(gene)

        # reads variants
        variants = self.queries.get_non_synonymous_variants_by_gene(gene_name)

        fig = dcc.Markdown("""**No variants for the current selection**""")
        if variants.shape[0] > 0:
            # reads total number of samples and calculates frequencies
            count_samples = self.queries.count_ena_samples()
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

            selected_variants_trace = None
            if selected_variants:
                selected_variants_trace = go.Scatter(
                    x=[int(v.get("dna_mutation").split(":")[0]) for v in selected_variants],
                    y=[v.get("frequency") for v in selected_variants],
                    name="selected variants",
                    mode='markers',
                    # opacity=0.5,
                    marker=dict(
                        symbol="circle",
                        color="blue",
                        size=10,
                        showscale=False
                    ),
                    xaxis='x',
                    showlegend=True,
                    text=["{} ({})".format(v.get("hgvs_p"), v.get("annotation")) for v in selected_variants],
                    hovertemplate='<b>%{text}</b><br>' +
                                  'Allele frequency: %{y:.5f}<br>' +
                                  'Genomic Position: %{x}'
                )

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
                xaxis='x',
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
                    xaxis='x',
                    legendgroup='domains'
                ))

            data = variants_traces + [gene_trace] + domain_traces
            if selected_variants_trace is not None:
                data.append(selected_variants_trace)
            layout = go.Layout(
                template="plotly_white",
                xaxis=dict(
                    domain=[0, 1.0],
                    tickformat=',d',
                    hoverformat=',d',
                    ticksuffix=" bp",
                    ticks="outside",
                    visible=True,
                    anchor="y2",
                    showspikes=True,
                    spikemode='across',
                    spikethickness=2
                ),
                xaxis2=dict(
                    domain=[0, 1.0],
                    anchor='y2',
                    visible=False
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
                margin=go.layout.Margin(l=0, r=0, b=0, t=20)
            )
            fig = go.Figure(data=data, layout=layout)

            # add vertical transparent rectangles with domains
            for d, c in zip(pfam_protein_features,
                            plotly.express.colors.qualitative.Plotly[0:len(pfam_protein_features)]):
                fig.add_vrect(
                    x0=gene_start + int(d["start"]),
                    x1=gene_start + int(d["end"]),
                    line_width=0, fillcolor=c, opacity=0.1)

            # add horizontal lines on the frequency boundaries
            fig.add_hline(y=0.1, line_width=1, line_dash="dash", line_color=VERY_COMMON_VARIANTS_COLOR)
            fig.add_hline(y=0.01, line_width=1, line_dash="dash", line_color=COMMON_VARIANTS_COLOR)
            fig.add_hline(y=0.001, line_width=1, line_dash="dash", line_color=RARE_VARIANTS_COLOR)

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
                    id="top-occurring-variants-table",
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
                    sort_by=[{"column_id": "frequency", "direction": "desc"}],
                    row_selectable='multi'
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

    def get_cooccurrence_heatmap(self, gene_name, selected_variants, metric="frequency", min_occurrences=5):
        data = self.queries.get_variants_cooccurrence_by_gene(gene_name=gene_name, min_cooccurrence=min_occurrences)
        fig = dcc.Markdown("""**No co-occurrent variants for the current selection**""")
        if data is not None and data.shape[0] > 0:

            # NOTE: transpose matrix manually as plotly transpose does not work with labels
            # the database return the upper diagonal, the lower is best for plots
            data.position_one = data.variant_one.transform(lambda x: int(x.split(":")[0]))
            data.reference_one = data.variant_one.transform(lambda x: x.split(":")[1].split(">")[0])
            data.alternate_one = data.variant_one.transform(lambda x: x.split(":")[1].split(">")[1])
            data.position_two = data.variant_two.transform(lambda x: int(x.split(":")[0]))
            data.reference_two = data.variant_two.transform(lambda x: x.split(":")[1].split(">")[0])
            data.alternate_two = data.variant_two.transform(lambda x: x.split(":")[1].split(">")[1])
            data.sort_values(["position_two", "reference_two", "alternate_two",
                              "position_one", "reference_two", "alternate_two"], inplace=True)

            # shortens very long ids
            def shorten_variant_id(variant_id):
                max_variant_length = 15
                if len(variant_id) <= max_variant_length:
                    return variant_id
                else:
                    return variant_id[0:max_variant_length] + "..."

            all_variants = [shorten_variant_id(v) for v in data.variant_one.dropna().unique()]

            values = np.array_split(data[metric], len(all_variants))
            texts = np.array_split(
                data[["hgvs_p_one", "hgvs_p_two"]].apply(
                    lambda x: "{} - {}".format(x.hgvs_p_one, x.hgvs_p_two), axis=1),
                len(all_variants))
            hovertemplate = '<b>%{text}</b><br>' + 'Cooccurrence: %{z:.5f}<br>' + 'Variant one: %{x}<br>' + 'Variant two: %{y}'
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
                values_selected = np.array_split(data[["variant_one", "variant_two", metric]].apply(
                    lambda x: x[metric] if x.variant_one in selected_variant_ids or
                                            x.variant_two in selected_variant_ids else None, axis=1),
                    len(all_variants))
                texts_selected = np.array_split(data[["variant_one", "variant_two", "hgvs_p_one", "hgvs_p_two"]].apply(
                    lambda x: "{} - {}".format(x.hgvs_p_one, x.hgvs_p_two)
                    if x.variant_one in selected_variant_ids or x.variant_two in selected_variant_ids else None,
                    axis=1),
                    len(all_variants))

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
                template="plotly_white",
                height=700,
                yaxis=dict(
                    visible=True,
                    tickfont={"size": 10},
                    showgrid=False,
                    showspikes=True,
                    spikemode='toaxis',
                    spikethickness=2
                ),
                xaxis=dict(
                    tickangle=-45,
                    tickfont={"size": 10},
                    showgrid=False,
                    showspikes=True,
                    spikemode='toaxis',
                    spikethickness=2
                ),
                margin=go.layout.Margin(l=0, r=0, b=0, t=0)
            )
            traces = [heatmap]
            if selected_variants:
                traces.append(heatmap_selected)
            fig = go.Figure(data=traces, layout=layout)

            # the y index is reversed in plotly heatmap
            fig.update_yaxes(autorange="reversed")

        return fig

    def get_variants_abundance_plot(self, bin_size=50, plotly_config=None):

        # reads genes and domains across the whole genome
        genes = sorted(self.queries.get_genes_metadata(), key=lambda x: int(x.data.get("start")))
        domains = []
        for g in genes:
            domains.extend([(g, d) for d in self.queries.get_pfam_domains(g)])

        # reads variants abundance
        variant_abundance = self.queries.get_variant_abundance_histogram(bin_size=bin_size)

        # reads conservation and bins it
        conservation = self.queries.get_conservation_table(bin_size=bin_size)

        # joins variant abundance and conservation
        data = variant_abundance.set_index("position_bin").join(conservation.set_index("position_bin"))
        data.reset_index(inplace=True)
        data.fillna(0, inplace=True)

        layout = go.Layout(
            template="plotly_white",
            xaxis=dict(
                domain=[0, 1.0],
                tickformat=',d',
                hoverformat=',d',
                ticksuffix=" bp",
                ticks="outside",
                visible=True,
                anchor="y7",
                showspikes=True,
                spikemode='across',
                spikethickness=2
            ),
            xaxis2=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            xaxis3=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            xaxis4=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            xaxis5=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            xaxis6=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            xaxis7=dict(
                domain=[0, 1.0],
                anchor='y7',
                visible=False
            ),
            yaxis=dict(
                domain=[0.9, 1.0],
                anchor='x7',
            ),
            yaxis2=dict(
                domain=[0.6, 0.9],
                anchor='x7'
            ),
            yaxis3=dict(
                domain=[0.45, 0.6],
                anchor='x7',
                visible=False,
            ),
            yaxis4=dict(
                domain=[0.3, 0.45],
                anchor='x7',
                visible=False,
            ),
            yaxis5=dict(
                domain=[0.15, 0.3],
                anchor='x7',
                visible=False,
            ),
            yaxis6=dict(
                domain=[0.05, 0.1],
                anchor='x7',
                visible=False,
            ),
            yaxis7=dict(
                domain=[0.0, 0.05],
                anchor='x7',
                visible=False,
            ),
            margin=go.layout.Margin(l=0, r=0, b=0, t=30),
            legend={'traceorder': 'normal'}
        )

        gene_traces = []
        for g, c in zip(genes, cycle(plotly.express.colors.sequential.Reds)):
            gene_start = int(g.data["start"])
            gene_end = int(g.data["end"])
            gene_name = g.data.get('name')
            gene_traces.append(go.Scatter(
                mode='lines',
                x=[gene_start, gene_end, gene_end, gene_start],
                y=[0, 0, 1, 1],
                name=gene_name,
                fill="toself",
                fillcolor=c,
                hovertext=gene_name,
                hoveron="fills",
                line=dict(width=0),
                yaxis='y6',
                xaxis='x',
                legendgroup='genes'
            ))

        domain_traces = []
        for (g, d), c in zip(domains, cycle(plotly.express.colors.sequential.Purples)):
            gene_start = int(g.data.get("start"))
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
                yaxis='y7',
                xaxis='x',
                legendgroup='domains',
                showlegend=False
            ))

        fig = go.Figure(
            data=[
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
                                line_color=plotly.express.colors.sequential.Blues[-4], line_width=1),
                     go.Scatter(x=data.position_bin, y=data.conservation, yaxis='y3',
                                text="Conservation SARS-CoV-2", textposition="top right", showlegend=False,
                                fill='tozeroy', line_color="grey", line_width=1),
                     go.Scatter(x=data.position_bin, y=data.conservation_sarbecovirus, yaxis='y4',
                                text="Conservation SARS-like betacoronavirus", textposition="top right",
                                showlegend=False, fill='tozeroy', line_color="grey", line_width=1),
                     go.Scatter(x=data.position_bin, y=data.conservation_vertebrates, yaxis='y5',
                                text="Conservation vertebrates", textposition="top right", showlegend=False,
                                fill='tozeroy', line_color="grey", line_width=1)
                 ] + gene_traces + domain_traces, layout=layout)

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
                config=plotly_config
            ),
            dcc.Markdown("""
                ***Genome view*** *representing the abundance of variants and ConsHMM (Arneson, 2019) conservation using a bin size of {} bp.*
                
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