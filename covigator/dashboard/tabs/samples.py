import functools

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Output, Input
from sqlalchemy.orm import Session
from covigator.dashboard.figures.samples import SampleFigures
from covigator.database.model import DataSource
from covigator.database.queries import Queries

DEFAULT_DATA_SOURCE = DataSource.ENA.name

ID_SLIDER_MIN_SAMPLES = 'slider-min-samples-per-country'
ID_DROPDOWN_DATA_SOURCE = "dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'dropdown-country'
ID_DROPDOWN_GENE = 'dropdown-gene-overall-mutations'
ID_ACCUMULATED_SAMPLES_GRAPH = 'accumulated-samples-per-country'
ID_DN_DS_GRAPH = 'dn_ds_graph'


@functools.lru_cache()
def get_tab_samples(queries: Queries):
    return dbc.Card(
        dbc.CardBody(
            children=[
                get_samples_tab_left_bar(queries),
                get_samples_tab_graphs()
        ])
    )


def get_samples_tab_graphs():
    return html.Div(
        className="ten columns",
        style={'overflow': 'scroll', "height": "900px"},
        children=[
            html.Br(),
            html.Div(id=ID_ACCUMULATED_SAMPLES_GRAPH),
            html.Br(),
            html.Div(id=ID_DN_DS_GRAPH),
        ])


def get_samples_tab_left_bar(queries: Queries):
    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""Select data source"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_DATA_SOURCE,
                options=[{'label': DataSource.ENA.name, 'value': DataSource.ENA.name},
                         {'label': DataSource.GISAID.name, 'value': DataSource.GISAID.name}],
                value=DEFAULT_DATA_SOURCE,
                clearable=False,
                multi=False
            ),
            html.Br(),
            dcc.Markdown("""Select one or more genes"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_GENE,
                options=[{'label': g.name, 'value': g.name} for g in queries.get_genes()],
                value=None,
                multi=True
            ),
            html.Br(),
            dcc.Markdown("""Select one or more countries"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_COUNTRY,
                options=[{'label': c, 'value': c} for c in queries.get_countries(DEFAULT_DATA_SOURCE)],
                value=None,
                multi=True
            ),
            html.Br(),
            dcc.Markdown("""**Accumulated samples by country**"""),
            dcc.Markdown("""Minimum number of samples per country"""),
            dcc.Slider(
                id=ID_SLIDER_MIN_SAMPLES,
                min=0,
                max=10000,
                step=100,
                value=100,
                dots=False,
                marks={i: '{}'.format(i) for i in [0, 2500, 5000, 7500, 10000]},
                tooltip=dict(always_visible=False, placement="right")
            ),
        ])


def set_callbacks_samples_tab(app, session: Session):

    queries = Queries(session=session)
    figures = SampleFigures(queries)

    @app.callback(
        Output(ID_DROPDOWN_COUNTRY, 'options'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'))
    def set_domains(source):
        """
        Updates the country drop down list when the data source is changed
        """
        return [{'label': c, 'value': c} for c in queries.get_countries(source)]

    @app.callback(
        Output(ID_ACCUMULATED_SAMPLES_GRAPH, 'children'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'),
        Input(ID_DROPDOWN_COUNTRY, 'value'),
        Input(ID_SLIDER_MIN_SAMPLES, 'value'),
        suppress_callback_exceptions=True
    )
    def update_accumulated_samples_by_country(data_source, countries, min_samples):
        return html.Div(children=figures.get_accumulated_samples_by_country_plot(
            data_source=data_source, countries=countries, min_samples=min_samples if countries is None else 0))

    @app.callback(
        Output(ID_DN_DS_GRAPH, 'children'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'),
        Input(ID_DROPDOWN_COUNTRY, 'value'),
        Input(ID_DROPDOWN_GENE, 'value'),
        suppress_callback_exceptions=True
    )
    def update_dn_ds_graph(data_source, countries, genes):
        return html.Div(children=figures.get_dnds_by_gene_plot(
            data_source=data_source, countries=countries, genes=genes))
