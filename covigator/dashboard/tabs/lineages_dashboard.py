from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from covigator.dashboard.figures import VARIANT_TYPE_COLOR_MAP
from covigator.dashboard.figures.figures import PLOTLY_CONFIG, MARGIN, TEMPLATE
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

def get_tab_lineages_dashboard():

    return dbc.CardBody(
            children=[
                get_lineages_dashboard_left_bar(),
                html.Div(
                    className="one column",
                    children=[html.Br()])
                #get_samples_tab_graphs()
        ])


def get_lineages_dashboard_left_bar():

    return html.Div(
        className="two columns",
        children=[
            html.Div([
                html.Div(dbc.Button(
                    "Lineage defining mutations",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'}
                )),
                html.Br(),
                html.Div(dbc.Button(
                    "Lineage annotation",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'})),
                html.Br(),
                html.Div(dbc.Button(
                    "Comparison of lineages",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'})),
                html.Br(),
                html.Div(dbc.Button(
                    "Geographical distribution",
                    color="secondary",
                    className="me-1",
                    style={'font-size': '100%'}))
            ]),
            html.Hr(),
            dcc.Markdown("""
                        Lineages assignments in the knowledge base are based on pangolin  nomenclature. Variants of concern were
                        further annotated with WHO label, variant of concern information and determining mutations as defined
                        in the cov-lineages [constellations](https://cov-lineages.org/constellations.html).
                         """, style={"font-size": 16})
        ])

def get_lineages_dashboard_graphs():
    pass