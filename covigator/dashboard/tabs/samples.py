import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
from covigator.dashboard.figures import Figures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_SLIDER_MIN_SAMPLES = 'slider-min-samples-per-country'
ID_DROPDOWN_DATA_SOURCE = "dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'dropdown-country'
ID_ACCUMULATED_SAMPLES_GRAPH = 'accumulated-samples-per-country'


def get_tab_samples(queries: Queries):
    return dcc.Tab(
        label="Samples",
        style=TAB_STYLE,
        selected_style=TAB_SELECTED_STYLE,
        children=[
            html.Div(
                id='ena-samples-body',
                className="row container-display",
                style={'overflow': 'scroll'}, # 'top': 0, 'bottom': 0, position: fixed
                children=[
                    get_samples_tab_left_bar(queries),
                    get_samples_tab_graphs()
                ])
        ]
    )


def get_samples_tab_graphs():
    return html.Div(
        className="ten columns",
        style={'overflow': 'scroll', "height": "900px"},
        children=[
            html.Br(),
            html.Div(id=ID_ACCUMULATED_SAMPLES_GRAPH),
            html.Br()
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
                value=None,
                multi=False
            ),
            html.Br(),
            dcc.Markdown("""Select a country"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_COUNTRY,
                options=[{'label': c, 'value': c} for c in queries.get_countries()],
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
                value=1000,
                dots=False,
                tooltip=dict(always_visible=True, placement="right")
            ),
            html.Br(),
        ])


def set_callbacks_samples_tab(app, figures: Figures, queries: Queries):
    @app.callback(
        Output(ID_ACCUMULATED_SAMPLES_GRAPH, 'children'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value'),
        Input(ID_DROPDOWN_COUNTRY, 'value'),
        Input(ID_SLIDER_MIN_SAMPLES, 'value'),
    )
    def update_accumulated_samples_by_country(data_source, countries, min_samples):
        return html.Div(children=figures.get_accumulated_samples_by_country_plot(
            data_source=data_source, countries=countries, min_samples=min_samples))
