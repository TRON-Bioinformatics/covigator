import functools

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session
from covigator.dashboard.figures.samples import SampleFigures
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_APPLY_BUTTOM = 's-apply-buttom'


ID_SLIDER_MIN_SAMPLES = 'slider-min-samples-per-country'
ID_DROPDOWN_DATA_SOURCE = "dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'dropdown-country'
ID_DROPDOWN_GENE = 'dropdown-gene-overall-mutations'
ID_ACCUMULATED_SAMPLES_GRAPH = 'accumulated-samples-per-country'
ID_DN_DS_GRAPH = 'dn_ds_graph'


@functools.lru_cache()
def get_tab_samples(queries: Queries, data_source: DataSource):
    return dbc.CardBody(
            children=[
                get_samples_tab_left_bar(queries, data_source),
                html.Div(
                    className="one column",
                    children=[html.Br()]),
                get_samples_tab_graphs()
        ])


def get_samples_tab_graphs():
    return html.Div(
        className="nine columns",
        children=[
            html.Br(),
            html.Div(id=ID_ACCUMULATED_SAMPLES_GRAPH),
            html.Br(),
            html.Div(id=ID_DN_DS_GRAPH),
        ])


def get_samples_tab_left_bar(queries: Queries, data_source: DataSource):
    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            html.Div(
                dcc.Dropdown(
                    id=ID_DROPDOWN_DATA_SOURCE,
                    options=[{'label': data_source.name, 'value': data_source.name}],
                    value=data_source.name,
                    clearable=False,
                    multi=False,
                    disabled=True
                ), style={'display': 'none'}),
            dcc.Markdown("""Select one or more countries"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_COUNTRY,
                options=[{'label': c, 'value': c} for c in queries.get_countries(data_source.name)],
                value=None,
                multi=True
            ),
            html.Br(),
            dcc.Markdown("""Minimum number of samples per country"""),
            dcc.Slider(
                id=ID_SLIDER_MIN_SAMPLES,
                min=0,
                max=10000,
                step=100,
                value=100,
                dots=False,
                marks={i: '{}'.format(i) for i in [0, 2500, 5000, 7500, 10000]},
                tooltip=dict(always_visible=True, placement="right")
            ),
            html.Br(),
            dcc.Markdown("""Select one or more genes to show the dN/dS ratio on protein domains"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_GENE,
                options=[{'label': g2, 'value': g2} for g2 in sorted([g.name for g in queries.get_genes()])],
                value=None,
                multi=True
            ),
            html.Br(),
            html.Button('Apply', id=ID_APPLY_BUTTOM),
        ])


def set_callbacks_samples_tab(app, session: Session):

    queries = Queries(session=session)
    figures = SampleFigures(queries)

    countries_ena = queries.get_countries(DataSource.ENA.name)
    countries_gisaid = queries.get_countries(DataSource.GISAID.name)

    @app.callback(
        Output(ID_DROPDOWN_COUNTRY, 'options'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'))
    def set_countries(source):
        """
        Updates the country drop down list when the data source is changed
        """
        countries = []
        if source == DataSource.ENA.name:
            countries = [{'label': c, 'value': c} for c in countries_ena]
        elif source == DataSource.GISAID.name:
            countries = [{'label': c, 'value': c} for c in countries_gisaid]
        return countries

    @app.callback(
        Output(ID_SLIDER_MIN_SAMPLES, 'disabled'),
        Input(ID_DROPDOWN_COUNTRY, 'value'))
    def disable_minimum_number_samples(countries):
        """
        Disables the minimum number of samples option when countries are provided
        """
        return countries is not None

    @app.callback(
        Output(ID_ACCUMULATED_SAMPLES_GRAPH, 'children'),
        [Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_SLIDER_MIN_SAMPLES, 'value')
        ],
        suppress_callback_exceptions=True
    )
    def update_accumulated_samples_by_country(_, data_source, countries, min_samples):
        return html.Div(children=figures.get_accumulated_samples_by_country_plot(
            data_source=data_source, countries=countries, min_samples=min_samples if countries is None else 0))

    @app.callback(
        Output(ID_DN_DS_GRAPH, 'children'),
        [Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_DROPDOWN_GENE, 'value'),
        ],
        suppress_callback_exceptions=True
    )
    def update_dn_ds_graph(_, data_source, countries, genes):
        return html.Div(children=figures.get_dnds_by_gene_plot(
            data_source=data_source, countries=countries, genes=genes))
