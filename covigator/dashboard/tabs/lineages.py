from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session
from covigator.dashboard.figures.lineages import LineageFigures
from covigator.dashboard.tabs import get_mini_container, print_number
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_APPLY_BUTTOM = 'lineages-apply-buttom'


ID_DROPDOWN_DATA_SOURCE = "lineages-dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'lineages-dropdown-country'
ID_DROPDOWN_LINEAGE = 'lineages-dropdown-lineage'
ID_LINEAGES_GRAPH = 'lineages-graph'
ID_LINEAGES_TABLE = 'lineages-table'


def get_tab_lineages(queries: Queries, data_source: DataSource):
    return dbc.CardBody(
            children=[
                get_lineages_tab_left_bar(queries, data_source),
                html.Div(
                    className="one column",
                    children=[html.Br()]),
                get_lineages_tab_graphs()
        ])


def get_lineages_tab_graphs():
    return html.Div(
        className="nine columns",
        children=[
            html.Br(),
            html.Div(id=ID_LINEAGES_GRAPH),
            html.Hr(),
            html.Br(),
            html.Div(id=ID_LINEAGES_TABLE),
        ])


def get_lineages_tab_left_bar(queries: Queries, data_source: DataSource):

    lineages = queries.get_lineages(source=data_source.name)

    return html.Div(
        className="two columns",
        children=[
            html.P("Lineage information is derived from the mutated sequence using Pangolin."),
            html.Br(),
            html.Div(
                html.Span(
                    children=[
                        get_mini_container(
                            title="Lineages",
                            value=print_number(len(lineages))
                        )
                        ])),
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
            dcc.Markdown("""Select one or more lineages"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_LINEAGE,
                options=[{'label': c, 'value': c} for c in queries.get_lineages(data_source.name)],
                value=None,
                multi=True
            ),
            html.Br(),
            html.P("Select a single lineage to explore its corresponding mutations."),
            html.Button('Apply', id=ID_APPLY_BUTTOM),
        ])


def set_callbacks_lineages_tab(app, session: Session):

    queries = Queries(session=session)
    figures = LineageFigures(queries)

    countries_ena = queries.get_countries(DataSource.ENA.name)
    countries_gisaid = queries.get_countries(DataSource.GISAID.name)
    lineages_ena = queries.get_lineages(DataSource.ENA.name)
    lineages_gisaid = queries.get_lineages(DataSource.GISAID.name)

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
        Output(ID_DROPDOWN_LINEAGE, 'options'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'))
    def set_lineages(source):
        """
        Updates the country drop down list when the data source is changed
        """
        lineages = []
        if source == DataSource.ENA.name:
            lineages = [{'label': c, 'value': c} for c in lineages_ena]
        elif source == DataSource.GISAID.name:
            lineages = [{'label': c, 'value': c} for c in lineages_gisaid]
        return lineages

    @app.callback(
        Output(ID_LINEAGES_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_DROPDOWN_LINEAGE, 'value'),
        ],
        suppress_callback_exceptions=True
    )
    def update_lineages_plot(_, data_source, countries, lineages):
        return html.Div(children=figures.get_lineages_plot(
            data_source=data_source,
            countries=countries,
            lineages=lineages))

    @app.callback(
        Output(ID_LINEAGES_TABLE, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_DROPDOWN_LINEAGE, 'value'),
        ],
        suppress_callback_exceptions=True
    )
    def update_lineages_table(_, data_source, countries, lineages):
        return html.Div(children=figures.get_lineages_variants_table(
            data_source=data_source, lineages=lineages, countries=countries))
