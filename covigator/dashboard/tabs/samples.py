import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input

from covigator.dashboard.figures import Figures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE
from covigator.database.queries import Queries


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
            html.Div(id='accumulated-samples-per-country'),
            html.Br()
        ])


def get_samples_tab_left_bar(queries: Queries):
    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""Select a country"""),
            dcc.Dropdown(
                id='dropdown-country',
                options=[{'label': c, 'value': c} for c in queries.get_ena_countries()],
                value=None,
                multi=True
            )
        ])


def set_callbacks_samples_tab(app, figures: Figures, queries: Queries):
    @app.callback(
        Output('accumulated-samples-per-country', 'children'),
        Input('dropdown-country', 'value'),
    )
    def update_accumulated_samples_by_country(countries):
        return html.Div(children=figures.get_accumulated_samples_by_country_plot(countries=countries))
