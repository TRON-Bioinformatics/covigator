from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
import plotly

from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
from covigator.dashboard.tabs import get_mini_container, print_number, print_date
from covigator.database.model import DataSource, SAMPLE_GISAID_TABLE_NAME, JobStatus
from covigator.database.queries import Queries
import pandas as pd
import plotly.express as px


def get_tab_dataset_gisaid(queries: Queries):
    count_samples = queries.count_samples(source=DataSource.GISAID.name)

    return dbc.CardBody(
        children=[
            get_gisaid_overview_tab_left_bar(queries, count_samples),
            html.Div(
                className="one column",
                children=[html.Br()]),
            get_gisaid_overview_tab_graphs(queries=queries, count_samples=count_samples)
        ])


def get_gisaid_overview_tab_graphs(queries, count_samples):
    return html.Div(
        className="nine columns",
        children=get_dataset_gisaid_tab_graphs(queries, count_samples=count_samples))


def get_gisaid_overview_tab_left_bar(queries: Queries, count_samples):
    count_variants = queries.count_variants(source=DataSource.GISAID.name)
    count_variant_observations = queries.count_variant_observations(source=DataSource.GISAID.name)
    date_of_first_gisaid_sample = queries.get_date_of_first_sample(source=DataSource.GISAID)
    date_of_most_recent_gisaid_sample = queries.get_date_of_most_recent_sample(source=DataSource.GISAID)

    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""
                The GISAID database (https://www.gisaid.org/) provides DNA assemblies, geographical information and
                other metadata about SARS-CoV-2 samples. 
                The processing pipeline runs alignment to the reference genome (bioypthon), 
                variant calling (custom code), normalization (vt and BCFtools), annotation (SnpEff) 
                and finally lineage determination (pangolin). 
                """, style={"font-size": 16}),
            html.Br(),
            html.Div(
                html.Span(
                    children=[
                        get_mini_container(
                            title="Samples",
                            value=print_number(count_samples)),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="First sample",
                            value=print_date(date_of_first_gisaid_sample)),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Latest sample",
                            value=print_date(date_of_most_recent_gisaid_sample)),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Unique mutations",
                            value=print_number(count_variants)),
                        html.Br(),
                        html.Br(),
                        get_mini_container(
                            title="Mutation calls",
                            value=print_number(count_variant_observations)),
                    ])
            ),
            html.Br(),
            dcc.Markdown("""
                    There are two sample exclusion criteria: 
                    * Horizontal coverage below 20 % of the reference genome
                    * A ratio of ambiguous bases greater than 0.2 
                    
                    Furthermore, all mutations involving ambiguous bases are excluded. 
                    """, style={"font-size": 16}),
        ])


def get_dataset_gisaid_tab_graphs(queries: Queries, count_samples):

    return html.Div(
        children=[
            html.Br(),
            html.Div(children=[
                html.Div(children=get_plot_coverage(queries),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
                html.Div(children=get_plot_bad_bases_ratio(queries, count_samples),
                         className="six columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
            ]),
        ])


def get_plot_coverage(queries: Queries):
    sql_query = """
    select count(*), (sequence_length::float / 29903 * 100)::int as coverage
    from {table} 
    where status = '{status}'
    group by coverage
    """.format(table=SAMPLE_GISAID_TABLE_NAME, status=JobStatus.FINISHED.name)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.bar(data_frame=data, x="coverage", y="count", log_y=True, color="count",
                 color_continuous_scale=plotly.colors.sequential.Brwnyl)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Num. of samples (log)"},
        xaxis={'title': "Horizontal coverage (%)"},
        legend={'title': None}
    )
    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Horizontal coverage (%)**
        
        Horizontal coverage is estimated as the sequence length divided by the reference genome length.
        Some samples have thus an horizontal coverage larger than 100 %. 
        """)
    ]


def get_plot_bad_bases_ratio(queries: Queries, count_samples):
    sql_query = """
    select count(*), ((count_n_bases + count_ambiguous_bases)::float / sequence_length * 100)::int as bad_bases_ratio
    from {table} 
    where status = '{status}'
    group by bad_bases_ratio
    """.format(table=SAMPLE_GISAID_TABLE_NAME, status=JobStatus.FINISHED.name)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.bar(data_frame=data, x="bad_bases_ratio", y="count", log_y=True, color="count",
                 color_continuous_scale=plotly.colors.sequential.Brwnyl)
    fig.update_layout(
        margin=MARGIN,
        template=TEMPLATE,
        yaxis={'title': "Num. of samples (log)"},
        xaxis={'title': "Ratio of N and ambiguous bases (%)"},
        legend={'title': None}
    )


    return [
        dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
        dcc.Markdown("""
        **Ratio of N and ambiguous bases (%)**

        The ratio of N and ambiguous bases over the whole sequence length measures the quality of the assembly sequence.
        {} % of samples have a ratio <= 5 %.
        """.format(round(float(data[data["bad_bases_ratio"] <= 5]["count"].sum()) / count_samples), 1))
    ]