from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
from covigator.dashboard.tabs import get_mini_container, print_number, print_date
from covigator.database.model import DataSource, SAMPLE_ENA_TABLE_NAME, VariantObservation, SampleEna, JobStatus
from covigator.database.queries import Queries
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go


LIBRARY_STRATEGY_COLOR_MAP = {
    "WGA": "#377eb8",
    "WGS": "#4daf4a",
    "RNA-Seq": "#984ea3",
    "Targeted-Capture": "#ff7f00",
}

LIBRARY_LAYOUT_COLOR_MAP = {
    "PAIRED": "#969696",
    "SINGLE": "#cccccc"
}


def get_tab_dataset_ena(queries: Queries):

    return dbc.CardBody(
        children=[
            get_ena_overview_tab_left_bar(queries),
            html.Div(
                className="one column",
                children=[html.Br()]),
            get_ena_overview_tab_graphs(queries=queries)
        ])


def get_ena_overview_tab_graphs(queries):
    return html.Div(
        className="nine columns",
        children=get_dataset_ena_tab_graphs(queries))


def get_ena_overview_tab_left_bar(queries: Queries):
    count_samples = queries.count_samples(source=DataSource.ENA.name)
    count_variants = queries.count_variants(source=DataSource.ENA.name)
    count_variant_observations = queries.count_variant_observations(source=DataSource.ENA.name)
    date_of_first_sample = queries.get_date_of_first_sample(source=DataSource.ENA)
    date_of_most_recent_sample = queries.get_date_of_most_recent_sample(source=DataSource.ENA)

    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""
                        The ENA database provides the raw reads and metadata from SARS-CoV-2 samples through the ENA 
                        Portal API (https://www.ebi.ac.uk/ena/portal/api/). 
                        FASTQ files containing the raw reads were downloaded and MD5 checked.
                        The processing pipeline runs reads trimming (fastp), alignment (BWA mem), 
                        coverage analysis (samtools), variant calling (LoFreq), normalization (vt and BCFtools),
                        annotation (SnpEff and VAFator) and finally lineage determination (pangolin). 
                         """, style={"font-size": 16}),
            html.Br(),
            html.Div(
                html.Span(
                    children=[
                        get_mini_container(
                            title="Samples",
                            value=print_number(count_samples)
                        ),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="First sample",
                            value=print_date(date_of_first_sample)
                        ),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Latest sample",
                            value=print_date(date_of_most_recent_sample)
                        ),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Unique mutations",
                            value=print_number(count_variants)
                        ),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Mutation calls",
                            value=print_number(count_variant_observations)
                        )
                    ]
                )
            ),
            html.Br(),
            dcc.Markdown("""
                There are four exclusion criteria:
                * Non-human host samples 
                * Non Illumina sequenced samples (ie: a large proportion of nanopore sequenced samples are excluded)
                * Horizontal coverage below 20 % of the reference genome
                * Mean mapping quality below 10 or mean base call quality below 10
                
                Mutations are classified according to their VAF into clonal mutations (0.8 <= VAF <= 1.0), 
                intrahost mutations (0.05 <= VAF < 0.8) and finally low frequency mutations (VAF < 0.05).
                """, style={"font-size": 16}),
        ])


def get_dataset_ena_tab_graphs(queries: Queries):
    return html.Div(
        children=[
            html.Br(),
            html.Div(children=[
                html.Div(children=get_plot_library_strategies(queries),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
                html.Div(children=get_plot_number_reads(queries),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
                html.Br(),
                html.Div(children=get_plot_coverage(queries),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
                html.Div(children=get_plot_qualities(queries),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
            ]),
        ])


def get_plot_library_strategies(queries: Queries):
    data = get_data_library_strategies(queries)
    fig = px.bar(data_frame=data, x="library_strategy", y="count", color="library_layout",
                 color_discrete_map=LIBRARY_LAYOUT_COLOR_MAP, barmode="group")
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Num. of samples"},
        xaxis={'title': None},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Library strategies**
        
        The number of samples for each library preparation strategy.
        Only Illumina samples are included, other technologies are excluded from CoVigator at the moment.
        """)
    ]


def get_data_library_strategies(queries):
    sql_query = """
    select count(*) as count, library_strategy, library_layout
    from {table} 
    where status = '{status}'
    group by library_strategy, library_layout""".format(
        table=SAMPLE_ENA_TABLE_NAME,
        status=JobStatus.FINISHED.name
    )
    data = pd.read_sql_query(sql_query, queries.session.bind)
    return data


def get_plot_number_reads(queries: Queries):
    data = get_data_number_reads(queries)

    fig = get_scatter_with_marginals(
        data=data,
        x="num_reads_after",
        y="num_reads_before",
        x_title="Num. of reads after trimming",
        y_title="Reported num. of reads")

    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Number of reads**
       
        The reported number of reads corresponds to the number reported in the ENA API.
        The number of reads after trimming corresponds to the number of reads after running fastp on the FASTQ files.
        No samples are excluded based on the number of reads.
        """)
    ]


def get_data_number_reads(queries):
    query = queries.session.query(
        SampleEna.run_accession, SampleEna.num_reads.label("num_reads_after"),
        SampleEna.read_count.label("num_reads_before"), SampleEna.library_strategy)\
        .filter(SampleEna.status == JobStatus.FINISHED.name)
    data = pd.read_sql_query(query.statement, queries.session.bind)
    return data


def get_scatter_with_marginals(data, x, y, x_title, y_title):
    fig = go.FigureWidget()
    for k, v in LIBRARY_STRATEGY_COLOR_MAP.items():
        fig.add_scattergl(
            x=data[data.library_strategy == k][x],
            y=data[data.library_strategy == k][y],
            mode='markers',
            marker=dict(color=v),
            name=k,
            opacity=0.5,
            hovertext=data[data.library_strategy == k].run_accession
        )
        fig.add_box(x=data[data.library_strategy == k][x], yaxis='y2', marker=dict(color=v),
                    showlegend=False)
        fig.add_box(y=data[data.library_strategy == k][y], xaxis='x2', marker=dict(color=v),
                    showlegend=False)
    fig.layout = dict(
        template=TEMPLATE,
        margin=MARGIN,
        yaxis={'title': y_title,
               'domain': [0, 0.85], 'showgrid': True, 'zeroline': False, 'type': 'log'},
        xaxis={'title': x_title,
               'domain': [0, 0.85], 'showgrid': True, 'zeroline': False, 'type': 'log'},
        legend={'title': None, 'itemsizing': 'constant'},
        bargap=0,
        xaxis2={'domain': [0.85, 1], 'showgrid': False, 'zeroline': False, 'visible': False, 'showticklabels': False},
        yaxis2={'domain': [0.85, 1], 'showgrid': False, 'zeroline': False, 'visible': False, 'showticklabels': False},
        hovermode='closest',
    )
    return fig


def get_plot_coverage(queries: Queries):
    data = get_data_coverage(queries)

    fig = get_scatter_with_marginals(
        data=data,
        x="mean_depth",
        y="coverage",
        x_title="Mean depth of coverage",
        y_title="Genome coverage (%)")

    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Coverage**

        The percentage of the genome covered by at least one read versus the mean depth of coverage.
        Samples covering less than 20 % of the genome are excluded.
        """)
    ]


def get_data_coverage(queries):
    query = queries.session.query(
        SampleEna.run_accession, SampleEna.library_strategy, SampleEna.coverage, SampleEna.mean_depth) \
        .filter(SampleEna.status == JobStatus.FINISHED.name)
    data = pd.read_sql_query(query.statement, queries.session.bind)
    return data


def get_plot_qualities(queries: Queries):
    data = get_data_qualities(queries)

    fig = get_scatter_with_marginals(
        data=data,
        x="mean_base_quality",
        y="mean_mapping_quality",
        x_title="Mean BCQ",
        y_title="Mean MQ")

    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Mean mapping quality (MQ) and Base Call Quality (BCQ)**
        
        Samples covering with either mapping quality < 10 or base call quality < 10 are excluded
        """)
    ]


def get_data_qualities(queries):
    query = queries.session.query(
        SampleEna.run_accession, SampleEna.library_strategy, SampleEna.mean_mapping_quality,
        SampleEna.mean_base_quality)\
        .filter(SampleEna.status == JobStatus.FINISHED.name)
    data = pd.read_sql_query(query.statement, queries.session.bind)
    return data


def get_plot_clonal_variants(queries: Queries):
    data = get_data_clonal_variants(queries)
    fig = px.scatter(data_frame=data, color="variant_type", y="vaf", x="dp",
                     log_x=True, log_y=True, opacity=0.5, hover_name="variant_id",
                     marginal_x='box', marginal_y='box',
                     color_discrete_map=VARIANT_TYPE_COLOR_MAP)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "VAF"},
        xaxis={'title': "DP"},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Variant Allele Frequency (VAF) and Depth of Coverage (DP)**

        Variants with VAF >= 0.8 are considered clonal 
        """)
    ]


def get_data_clonal_variants(queries):
    query = queries.session.query(
        VariantObservation.variant_id, VariantObservation.vaf, VariantObservation.dp, VariantObservation.variant_type)
    data = pd.read_sql_query(query.statement, queries.session.bind)
    return data


def get_dataset_ena_tab_left_bar(queries: Queries):
    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""Select one or more genes"""),
            dcc.Dropdown(
                id="sample-dropdown",
                options=[{'label': g.name, 'value': g.name} for g in queries.get_genes()],
                value=None,
                multi=True
            )
        ])
