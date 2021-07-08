import dash_core_components as dcc
import dash_html_components as html
import plotly

from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE, get_mini_container
from covigator.database.model import DataSource, VariantType, SAMPLE_GISAID_TABLE_NAME
from covigator.database.queries import Queries
import pandas as pd
import plotly.express as px


def get_tab_dataset_gisaid(queries: Queries):

    count_samples = queries.count_samples(source=DataSource.GISAID.name)
    count_variants = queries.count_variant_observations(source=DataSource.GISAID.name)
    date_of_first_gisaid_sample = queries.get_date_of_first_sample(source=DataSource.GISAID)
    date_of_most_recent_gisaid_sample = queries.get_date_of_most_recent_sample(source=DataSource.GISAID)

    return dcc.Tab(
        label="GISAID dataset",
        style=TAB_STYLE,
        selected_style=TAB_SELECTED_STYLE,
        children=[
            html.Div(id="something-gisaid", children=[
                html.Div(className="two columns", children=[html.Br()]),
                html.Div(className="eight columns", children=[
                    html.Br(),
                    dcc.Markdown("""
                        The GISAID dataset was manually downloaded from the site https://www.gisaid.org/.
                        DNA assemblies and metadata were matched together and variant calling was done after performing 
                        a global alignment to the reference genome.
                        """),
                    html.Div(
                        className="row container-display",
                        children=[
                            get_mini_container(title="No. of samples", value=count_samples),
                            get_mini_container(title="No. of variant calls", value=count_variants),
                            get_mini_container(title="First sample", value=date_of_first_gisaid_sample),
                            get_mini_container(title="Latest sample", value=date_of_most_recent_gisaid_sample),
                        ]
                    )
                ]),
                html.Div(className="two columns", children=[html.Br()]),
            ]),
            html.Br(),
            html.Div(
                id='gisaid-dataset-body',
                className="row container-display twelve columns",
                children=[
                    get_dataset_gisaid_tab_graphs(queries=queries, count_samples=count_samples)
                ])
        ]
    )


def get_dataset_gisaid_tab_graphs(queries: Queries, count_samples):

    return html.Div(
        className="twelve columns",
        children=[
            html.Br(),
            html.Div(children=[
                html.Div(className="two columns", children=[html.Br()]),
                html.Div(className="four columns", children=get_plot_coverage(queries)),
                html.Div(className="four columns", children=get_plot_bad_bases_ratio(queries, count_samples)),
            ]),
        ])


def get_plot_coverage(queries: Queries):
    sql_query = """
    select count(*), (sequence_length::float / 29903 * 100)::int as coverage
    from {table} 
    where finished
    group by coverage
    """.format(table=SAMPLE_GISAID_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.bar(data_frame=data, x="coverage", y="count", log_y=True, color="count",
                 color_continuous_scale=plotly.colors.sequential.Brwnyl)
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
        **Horizontal coverage (%)**
        
        Horizontal coverage is estimated as the sequence length divided by the reference genome length.
        Some samples have thus an horizontal coverage larger than 100 %.
        Samples covering less than 20% of the whole genome are excluded. 
        """)
    ]


def get_plot_bad_bases_ratio(queries: Queries, count_samples):
    sql_query = """
    select count(*), ((count_n_bases + count_ambiguous_bases)::float / sequence_length * 100)::int as bad_bases_ratio
    from {table} 
    where finished
    group by bad_bases_ratio
    """.format(table=SAMPLE_GISAID_TABLE_NAME)
    data = pd.read_sql_query(sql_query, queries.session.bind)
    fig = px.bar(data_frame=data, x="bad_bases_ratio", y="count", log_y=True, color="count",
                 color_continuous_scale=plotly.colors.sequential.Brwnyl)
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
        **Ratio of N and ambiguous bases (%)**

        The ratio of N and ambiguous bases over the whole sequence length measures the quality of the assembly sequence.
        {} % of samples have a ratio <= 5 %.
        Samples with a ratio higher than 0.2 are excluded.
        All variant calls that contain a N or an ambiguous base are filtered out.
        """.format(round(float(data[data["bad_bases_ratio"] <= 5]["count"].sum()) / count_samples), 1))
    ]