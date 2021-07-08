import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
from sqlalchemy.orm import Session
from covigator.dashboard.figures.samples import SampleFigures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE
from covigator.database.model import DataSource, VariantType
from covigator.database.queries import Queries


def get_tab_dataset_gisaid(queries: Queries):
    return dcc.Tab(
        label="GISAID dataset",
        style=TAB_STYLE,
        selected_style=TAB_SELECTED_STYLE,
        children=[
            html.Div(
                id='gisaid-dataset-body',
                className="row container-display",
                style={'overflow': 'scroll'}, # 'top': 0, 'bottom': 0, position: fixed
                children=[
                    get_dataset_gisaid_tab_left_bar(queries),
                    get_dataset_gisaid_tab_graphs()
                ])
        ]
    )


def get_dataset_gisaid_tab_graphs():
    return html.Div(
        className="ten columns",
        style={'overflow': 'scroll', "height": "900px"},
        children=[
            html.Br(),
            html.Div(id="something-gisaid"),
            html.Br(),
            html.Div(children=[
                html.Div(id="else-gisaid", className="five columns"),
                html.Div(id="another-gisaid", className="five columns"),
            ]),
        ])


def get_dataset_gisaid_tab_left_bar(queries: Queries):
    return html.Div(
        className="two columns",
        children=[
            html.Br(),
            dcc.Markdown("""Select one or more genes"""),
            dcc.Dropdown(
                id="sample-dropdown-gisaid",
                options=[{'label': g.name, 'value': g.name} for g in queries.get_genes()],
                value=None,
                multi=True
            )
        ])


def set_callbacks_dataset_gisaid_tab(app, session: Session):

    queries = Queries(session=session)
    figures = SampleFigures(queries)

    #@app.callback(
    #    Output(ID_ACCUMULATED_SAMPLES_GRAPH, 'children'),
    #    Input(ID_DROPDOWN_DATA_SOURCE, 'value'),
    #    Input(ID_DROPDOWN_COUNTRY, 'value'),
    #    Input(ID_SLIDER_MIN_SAMPLES, 'value'),
    #)
    #def update_accumulated_samples_by_country(data_source, countries, min_samples):
    #    return html.Div(children=figures.get_accumulated_samples_by_country_plot(
    #        data_source=data_source, countries=countries, min_samples=min_samples))
