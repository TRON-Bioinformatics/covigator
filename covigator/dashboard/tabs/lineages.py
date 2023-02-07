from datetime import timedelta, datetime
from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session
from covigator.dashboard.figures.lineages import LineageFigures
from covigator.dashboard.tabs import get_mini_container, print_number, MONTH_PATTERN
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_APPLY_BUTTOM = 'lineages-apply-buttom'


ID_DROPDOWN_DATA_SOURCE = "lineages-dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'lineages-dropdown-country'
ID_DROPDOWN_LINEAGE = 'lineages-dropdown-lineage'
ID_DROPDOWN_PERIOD = 'lineages-dropdown-period'
ID_DROPDOWN_LINEAGE_DATE_RANGE_END = 'lineages-dropdown-date-range-end'
ID_DROPDOWN_LINEAGE_DATE_RANGE_START = 'lineages-dropdown-date-range-start'
ID_DROPDOWN_LINEAGE_DATE_RANGE_END_DIV = 'lineages-dropdown-date-range-end-div'
ID_SLIDER_PREVALENCE = 'lineages-slider-prevalence'
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
    # Get all available months from collection_date column in sample table
    months = queries.get_sample_months(MONTH_PATTERN, data_source=data_source.name)
    today = datetime.now()
    today_formatted = today.strftime(MONTH_PATTERN)
    oneyearago = today - timedelta(days=356)
    oneyearago_formatted = oneyearago.strftime(MONTH_PATTERN)

    return html.Div(
        className="two columns",
        children=[
            html.Div([
                html.Div(dbc.Button(
                    "Lineages accumulation",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'}
                )),
                html.Br(),
                html.Div(dbc.Button(
                    "Instantaneous lineage prevalence",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'}
                ))]),
            html.Hr(),
            html.P("Lineage information is derived from the consensus sequence using Pangolin."),
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
            dcc.Markdown("""Select time period for smoothing"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_PERIOD,
                options=[{'label': '{} days'.format(c), 'value': c} for c in range(1, 32)] + [{'label': 'Disable', 'value': False}],
                value=14,
                multi=False
            ),
            html.Br(),
            dcc.Markdown("""Select a start and end date"""),
            html.Div(children=[
                dcc.Dropdown(
                    id=ID_DROPDOWN_LINEAGE_DATE_RANGE_START,
                    options=[{'label': c, 'value': c} for c in months],
                    value=oneyearago_formatted,
                    multi=False,
                    clearable=False
                ),
                html.Div(
                    id=ID_DROPDOWN_LINEAGE_DATE_RANGE_END_DIV,
                    children=dcc.Dropdown(
                        id=ID_DROPDOWN_LINEAGE_DATE_RANGE_END,
                        options=[{'label': c, 'value': c} for c in months],
                        value=today_formatted,
                        multi=False,
                        clearable=False
                    ))]),
            html.Br(),
            dcc.Markdown("""Minimum prevalence of lineage in the time interval to be plotted"""),
            dcc.Slider(
                id=ID_SLIDER_PREVALENCE,
                min=0.0,
                max=0.2,
                step=0.01,
                value=0.01,
                marks={i: '{}'.format(i) for i in [0, 0.05, 0.1, 0.15, 0.2]},
                tooltip=dict(always_visible=False, placement="right")
            ),
            html.Br(),
            html.P("Select a single lineage to explore its corresponding mutations."),
            html.Button('Apply', id=ID_APPLY_BUTTOM),
        ])


def set_callbacks_lineages_tab(app, session: Session):

    queries = Queries(session=session)
    figures = LineageFigures(queries)

    countries_ena = queries.get_countries(DataSource.ENA.name)
    countries_covid19_portal = queries.get_countries(DataSource.COVID19_PORTAL.name)
    lineages_ena = queries.get_lineages(DataSource.ENA.name)
    lineages_covid19_portal = queries.get_lineages(DataSource.COVID19_PORTAL.name)

    # Get months from ENA/Covid19 table
    months_from_db_ena = queries.get_sample_months(MONTH_PATTERN, data_source=DataSource.ENA.name)
    months_from_db_covid19_portal = queries.get_sample_months(MONTH_PATTERN, data_source=DataSource.COVID19_PORTAL.name)

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
        elif source == DataSource.COVID19_PORTAL.name:
            countries = [{'label': c, 'value': c} for c in countries_covid19_portal]
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
        elif source == DataSource.COVID19_PORTAL.name:
            lineages = [{'label': c, 'value': c} for c in lineages_covid19_portal]
        return lineages

    @app.callback(
        Output(ID_DROPDOWN_LINEAGE_DATE_RANGE_END_DIV, 'children'),
        Input(ID_DROPDOWN_LINEAGE_DATE_RANGE_START, 'value'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value')
    )
    def update_dropdown_end_date(start_date, data_source):
        today = datetime.now()
        today_formatted = today.strftime(MONTH_PATTERN)
        months = []
        if data_source == DataSource.ENA.name:
            months = [m for m in months_from_db_ena if m >= start_date]
        elif data_source == DataSource.COVID19_PORTAL.name:
            months = [m for m in months_from_db_covid19_portal if m >= start_date]
        return dcc.Dropdown(
            id=ID_DROPDOWN_LINEAGE_DATE_RANGE_END,
            options=[{'label': c, 'value': c} for c in months],
            value=today_formatted,
            multi=False,
            clearable=False
        )

    @app.callback(
        Output(ID_LINEAGES_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_DROPDOWN_LINEAGE, 'value'),
            State(ID_DROPDOWN_PERIOD, 'value'),
            State(ID_DROPDOWN_LINEAGE_DATE_RANGE_START, 'value'),
            State(ID_DROPDOWN_LINEAGE_DATE_RANGE_END, 'value'),
            State(ID_SLIDER_PREVALENCE, 'value')
        ],
        suppress_callback_exceptions=True)
    def update_lineages_plot(_, data_source, countries, lineages, time_period, date_start, date_end, prevalence):
        return html.Div(children=figures.get_lineages_plot(
            data_source=data_source,
            countries=countries,
            lineages=lineages,
            time_period=time_period,
            date_start=date_start,
            date_end=date_end,
            prevalence=prevalence
        ))

    @app.callback(
        Output(ID_LINEAGES_TABLE, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_COUNTRY, 'value'),
            State(ID_DROPDOWN_LINEAGE, 'value'),
        ],
        suppress_callback_exceptions=True)
    def update_lineages_table(_, data_source, countries, lineages):
        return html.Div(children=figures.get_lineages_variants_table(
            data_source=data_source, lineages=lineages, countries=countries))
