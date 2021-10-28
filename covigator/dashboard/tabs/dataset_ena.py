import functools
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
from covigator.dashboard.tabs import get_mini_container, print_number, print_date
from covigator.database.model import DataSource, SAMPLE_ENA_TABLE_NAME, JOB_ENA_TABLE_NAME, \
    VARIANT_OBSERVATION_TABLE_NAME
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


@functools.lru_cache()
def get_tab_dataset_ena(queries: Queries):

    count_samples = queries.count_samples(source=DataSource.ENA.name)
    count_variants = queries.count_variant_observations(source=DataSource.ENA.name)
    date_of_first_sample = queries.get_date_of_first_sample(source=DataSource.ENA)
    date_of_most_recent_sample = queries.get_date_of_most_recent_sample(source=DataSource.ENA)

    return dbc.CardBody(
            children=[
                dcc.Markdown("""
                The ENA dataset and its metadata was downloaded using the ENA Portal API 
                (https://www.ebi.ac.uk/ena/portal/api/). FASTQ files containing the raw reads were downloaded from the
                 provided URLs for each sample. FASTQ files were MD5 checked after download. 
                 All non-human host samples were excluded.
                 """, style={"font-size": 16}),
                html.Br(),
                html.Div(
                    html.Span(
                        children=[
                            get_mini_container(
                                title="Samples",
                                value=print_number(count_samples)
                            ),
                            get_mini_container(
                                title="Variant calls",
                                value=print_number(count_variants)
                            ),
                            get_mini_container(
                                title="First sample",
                                value=print_date(date_of_first_sample)
                            ),
                            get_mini_container(
                                title="Latest sample",
                                value=print_date(date_of_most_recent_sample)
                            )
                        ]
                    )
                ),
                get_dataset_ena_tab_graphs(queries)
            ]
        )


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
    where finished
    group by library_strategy, library_layout""".format(table=SAMPLE_ENA_TABLE_NAME)
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
    sql_query = """
    select s.run_accession as run_accession, j.num_reads as num_reads_after, s.read_count as num_reads_before, s.library_strategy
    from {samples_table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where s.finished""".format(samples_table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
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
    sql_query = """
    select s.run_accession, s.library_strategy, j.coverage, j.mean_depth
    from {table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where finished""".format(table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
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
    sql_query = """
    select s.run_accession, s.library_strategy, j.mean_mapping_quality, j.mean_base_quality
    from {table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where finished""".format(table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
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
    sql_query = """
    select variant_id, vaf, dp, variant_type
    from {table} 
    where source='ENA'""".format(table=VARIANT_OBSERVATION_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
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
