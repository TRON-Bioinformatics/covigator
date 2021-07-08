import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
from sqlalchemy.orm import Session

from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE, get_mini_container, COLOR_OVERVIEW_MINI_CONTAINER, \
    print_number, print_date
from covigator.database.model import DataSource, VariantType, SAMPLE_ENA_TABLE_NAME, JOB_ENA_TABLE_NAME, \
    VARIANT_OBSERVATION_TABLE_NAME
from covigator.database.queries import Queries
import pandas as pd
import plotly.express as px


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

    count_samples = queries.count_samples(source=DataSource.ENA.name)
    count_variants = queries.count_variant_observations(source=DataSource.ENA.name)
    date_of_first_sample = queries.get_date_of_first_sample(source=DataSource.ENA)
    date_of_most_recent_sample = queries.get_date_of_most_recent_sample(source=DataSource.ENA)

    return dcc.Tab(
        label="ENA dataset",
        style=TAB_STYLE,
        selected_style=TAB_SELECTED_STYLE,
        children=[
            html.Div(id="something-ena", children=[
                html.Div(className="one columns", children=[html.Br()]),
                html.Div(className="nine columns", children=[
                    html.Br(),
                    dcc.Markdown("""
                                The ENA dataset was downloaded using the API https://www.ebi.ac.uk/ena/portal/api/.
                                Metadata was downloaded from the API and the FASTQ files with the raw reads were 
                                downloaded from the provided URLs for each sample. 
                                FASTQ files were MD5 checked after download.
                                All samples were the host was not human were excluded.
                                """, style={"font-size": 16}),
                    get_dataset_ena_tab_graphs(queries)
                ]),
                html.Div(html.Br(), className="one columns"),
                html.Div(className="one columns", children=[
                    html.Br(),
                    get_mini_container(
                        title="Samples",
                        value=print_number(count_samples),
                        color=COLOR_OVERVIEW_MINI_CONTAINER
                    ),
                    get_mini_container(
                        title="Variant calls",
                        value=print_number(count_variants),
                        color=COLOR_OVERVIEW_MINI_CONTAINER
                    ),
                    get_mini_container(
                        title="First sample",
                        value=print_date(date_of_first_sample),
                        color=COLOR_OVERVIEW_MINI_CONTAINER
                    ),
                    get_mini_container(
                        title="Latest sample",
                        value=print_date(date_of_most_recent_sample),
                        color=COLOR_OVERVIEW_MINI_CONTAINER
                    )
                ]),
            ]),
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
    sql_query = """
    select count(*) as count, library_strategy, library_layout
    from {table} 
    where finished
    group by library_strategy, library_layout""".format(table=SAMPLE_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.bar(data_frame=data, x="library_strategy", y="count", color="library_layout",
                 color_discrete_map=LIBRARY_LAYOUT_COLOR_MAP)
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


def get_plot_number_reads(queries: Queries):
    sql_query = """
    select s.run_accession as run_accession, j.num_reads as num_reads_after, s.read_count as num_reads_before, s.library_strategy
    from {samples_table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where s.finished""".format(samples_table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.scatter(data_frame=data, x="num_reads_after", y="num_reads_before", color="library_strategy",
                     opacity=0.5, hover_name="run_accession", log_y=True, log_x=True,
                     marginal_x='box', marginal_y='box',
                     color_discrete_map=LIBRARY_STRATEGY_COLOR_MAP)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Reported num. of reads"},
        xaxis={'title': "Num. of reads after trimming"},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Number of reads**
       
        The reported number of reads corresponds to the number reported in the ENA API.
        The number of reads after trimming corresponds to the number of reads after running fastp on the FASTQ files.
        No samples are excluded based on the number of reads.
        """)
    ]


def get_plot_coverage(queries: Queries):
    sql_query = """
    select s.run_accession, s.library_strategy, j.coverage, j.mean_depth
    from {table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where finished""".format(table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.scatter(data_frame=data, color="library_strategy", y="coverage", x="mean_depth",
                     log_x=True, log_y=True, opacity=0.5, hover_name="run_accession",
                     marginal_x='box', marginal_y='box',
                     color_discrete_map=LIBRARY_STRATEGY_COLOR_MAP)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Genome coverage (%)"},
        xaxis={'title': "Mean depth of coverage"},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Horizontal coverage (%)**

        The percentage of the genome covered by at least one read.
        Samples covering less than 20 % of the genome are excluded.
        """)
    ]


def get_plot_qualities(queries: Queries):
    sql_query = """
    select s.run_accession, s.library_strategy, j.mean_mapping_quality, j.mean_base_quality
    from {table} as s join {jobs_table} as j on s.run_accession = j.run_accession 
    where finished""".format(table=SAMPLE_ENA_TABLE_NAME, jobs_table=JOB_ENA_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    #fig_mq = px.box(data_frame=data, color="library_strategy", y="mean_mapping_quality",
    #                color_discrete_map=LIBRARY_STRATEGY_COLOR_MAP)
    #fig_bcq = px.box(data_frame=data, color="library_strategy", y="mean_base_quality",
    #                   color_discrete_map=LIBRARY_STRATEGY_COLOR_MAP)
    fig = px.scatter(data_frame=data, color="library_strategy", y="mean_mapping_quality", x="mean_base_quality",
                     log_x=True, log_y=True, opacity=0.5, hover_name="run_accession",
                     marginal_x='box', marginal_y='box',
                     color_discrete_map=LIBRARY_STRATEGY_COLOR_MAP)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Mean MQ"},
        xaxis={'title': "Mean BCQ"},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Mean mapping quality and Base call quality**
        
        Samples covering with either mapping quality < 10 or base call quality < 10 are excluded
        """)
    ]


def get_plot_clonal_variants(queries: Queries):
    sql_query = """
    select variant_id, vaf, dp, variant_type
    from {table} 
    where source='ENA'""".format(table=VARIANT_OBSERVATION_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
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
