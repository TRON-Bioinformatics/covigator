# Run this app with `python dashboard.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from dash.dependencies import Output, Input

from covigator.model import EnaRun
from covigator.database import Database

import flask
import pkg_resources
import os
import random


# creates the Dash application
app = dash.Dash(
    __name__,
    title="CoVigator",
    meta_tags=[
        # A description of the app, used by e.g.
        # search engines when displaying search results.
        {
            'name': 'description',
            'content': 'CoVigator - monitoring Sars-Cov-2 mutations'
        },
        # A tag that tells Internet Explorer (IE)
        # to use the latest renderer version available
        # to that browser (e.g. Edge)
        {
            'http-equiv': 'X-UA-Compatible',
            'content': 'IE=edge'
        },
        # A tag that tells the browser not to scale
        # desktop widths to fit mobile screens.
        # Sets the width of the viewport (browser)
        # to the width of the device, and the zoom level
        # (initial scale) to 1.
        #
        # Necessary for "true" mobile support.
        {
          'name': 'viewport',
          'content': 'width=device-width, initial-scale=1.0'
        }
    ],
    #external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css']
)

# sets up page layout
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

logo = "/assets/CoVigator_logo_txt.png"

# initialize the database connection
session = Database().get_database_session()

# read all samples data
runs = pd.read_sql(session.query(EnaRun).statement, session.bind)


samples_by_date_and_country = runs[["first_created", "run_accession", "country_alpha_3"]].groupby(["first_created", "country_alpha_3"]).count()
samples_by_date_and_country.reset_index(inplace=True)
samples_by_date_and_country.rename(columns={"run_accession": "count"}, inplace=True)
samples_by_date_and_country['cumsum'] = samples_by_date_and_country.groupby(['country_alpha_3'])['count'].cumsum()
other_countries = samples_by_date_and_country[samples_by_date_and_country["cumsum"] < 30].reset_index().country_alpha_3
#dates = samples_by_date_and_country.first_created.unique()
#countries = samples_by_date_and_country.sort_values("cumsum", ascending=False).country_alpha_3.unique()
#index = pd.MultiIndex.from_product([dates, countries], names = ["first_created", "country_alpha_3"])
#data = pd.DataFrame(index = index).reset_index()
#data["count"] = 0
#data.set_index(["first_created", "country_alpha_3"]).join(samples_by_date_and_country, on=["first_created", "country_alpha_3"])
#data2 = data.set_index(["first_created", "country_alpha_3"]) + samples_by_date_and_country.set_index(["first_created", "country_alpha_3"])
#data2.fillna(0, inplace=True)
#data2.reset_index(inplace=True)
#data2['cumsum'] = data2.groupby(['country_alpha_3'])['count'].cumsum()
#data2.rename(columns={"first_created": "date", "country_alpha_3": "country"}, inplace=True)
#data3 = data2.groupby("country").max()
#other_countries = data3[data3["cumsum"] < 30].reset_index().country

# merge countries with less than 30 samples into OTHER
runs["country_merged"] = runs.country_alpha_3.transform(lambda c: "Other" if c in list(other_countries) or c is None or c == "None" else c)
samples_by_date_and_country_merged = runs[["first_created", "run_accession", "country_merged"]].groupby(["first_created", "country_merged"]).count()
samples_by_date_and_country_merged.reset_index(inplace=True)
samples_by_date_and_country_merged.rename(columns={"run_accession": "count"}, inplace=True)
samples_by_date_and_country_merged['cumsum'] = samples_by_date_and_country_merged.groupby(['country_merged'])['count'].cumsum()

# creates empty table with all dates and all countries
dates = samples_by_date_and_country_merged.first_created.unique()
countries = samples_by_date_and_country_merged.sort_values("cumsum", ascending=False).country_merged.unique()
index = pd.MultiIndex.from_product([dates, countries], names = ["first_created", "country_merged"])
data_merged = pd.DataFrame(index = index).reset_index()
data_merged["count"] = 0

# adds values into empty table
data2_merged = data_merged.set_index(["first_created", "country_merged"]) + samples_by_date_and_country_merged.set_index(["first_created", "country_merged"])
data2_merged.fillna(0, inplace=True)
data2_merged.reset_index(inplace=True)
data2_merged['cumsum'] = data2_merged.groupby(['country_merged'])['count'].cumsum()
data2_merged.rename(columns={"first_created": "date", "country_merged": "country"}, inplace=True)

# Create app layout
app.layout = html.Div(
    [
        dcc.Store(id="aggregate_data"),
        # empty Div to trigger javascript file for graph resizing
        html.Div(id="output-clientside"),
        html.Div(
            [
                html.Div(
                    [
                        html.Img(
                            src=logo,
                            id="covigator-logo",
                            style={
                                "height": "60px",
                                "width": "auto",
                                "margin-bottom": "25px",
                            },
                        )
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H3(
                                    "Monitoring Sars-Cov-2 mutations",
                                    style={"margin-bottom": "0px"},
                                ),
                                html.P("Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, "
                                       "necessitating preventive or therapeutic strategies and first steps towards an end to "
                                       "this pandemic were done with the approval of the first mRNA vaccines against SARS-CoV-2. "
                                       "We want to provide an interactive view on different types of variants that can be "
                                       "incorporated in global efforts to sustainably prevent or treat infections. Thus, we "
                                       "envision to help guiding global vaccine design efforts to overcome the threats of this "
                                       "pandemic.",
                                       style={"text-align": "left"},),
                                html.A(
                                    html.Button(
                                        "https://doi.org/10.1101/2021.02.04.429765",
                                        id="tron-button"),
                                        href="https://doi.org/10.1101/2021.02.04.429765",
                                )
                            ]
                        )
                    ],
                    #className="one-half column",
                    className="two column",
                    id="title",
                ),
                # html.Div(
                #     [
                #         html.P("Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, "
                #                "necessitating preventive or therapeutic strategies and first steps towards an end to "
                #                "this pandemic were done with the approval of the first mRNA vaccines against SARS-CoV-2. "
                #                "We want to provide an interactive view on different types of variants that can be "
                #                "incorporated in global efforts to sustainably prevent or treat infections. Thus, we "
                #                "envision to help guiding global vaccine design efforts to overcome the threats of this "
                #                "pandemic."),
                #         html.A(
                #             html.Button(
                #                 "https://doi.org/10.1101/2021.02.04.429765",
                #                 id="tron-button"),
                #             href="https://doi.org/10.1101/2021.02.04.429765",
                #         )
                #     ],
                #     className="one-third column",
                #     id="button",
                # ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.P(
                            "Filter by construction date (or select range in histogram):",
                            className="control_label",
                        ),
                        dcc.RangeSlider(
                            id="year_slider",
                            min=1960,
                            max=2017,
                            value=[1990, 2010],
                            className="dcc_control",
                        ),
                        html.P("Filter by well status:", className="control_label"),
                        dcc.RadioItems(
                            id="well_status_selector",
                            options=[
                                {"label": "All ", "value": "all"},
                                {"label": "Active only ", "value": "active"},
                                {"label": "Customize ", "value": "custom"},
                            ],
                            value="active",
                            labelStyle={"display": "inline-block"},
                            className="dcc_control",
                        ),
                        dcc.Dropdown(
                            id="well_statuses",
                            options=[],
                            multi=True,
                            value=[],
                            className="dcc_control",
                        ),
                        dcc.Checklist(
                            id="lock_selector",
                            options=[{"label": "Lock camera", "value": "locked"}],
                            className="dcc_control",
                            value=[],
                        ),
                        html.P("Filter by well type:", className="control_label"),
                        dcc.RadioItems(
                            id="well_type_selector",
                            options=[
                                {"label": "All ", "value": "all"},
                                {"label": "Productive only ", "value": "productive"},
                                {"label": "Customize ", "value": "custom"},
                            ],
                            value="productive",
                            labelStyle={"display": "inline-block"},
                            className="dcc_control",
                        ),
                        dcc.Dropdown(
                            id="well_types",
                            options=[],
                            multi=True,
                            value=[],
                            className="dcc_control",
                        ),
                    ],
                    className="pretty_container four columns",
                    id="cross-filter-options",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    [html.H6(runs.shape[0]), html.P("No. of Samples")],
                                    className="mini_container",
                                ),
                                html.Div(
                                    [html.H6(len(runs.country.unique())), html.P("No. of countries")],
                                    className="mini_container",
                                ),
                                html.Div(
                                    [html.H6(12345), html.P("No. of unique variants")],
                                    className="mini_container",
                                ),
                                html.Div(
                                    [html.H6(12345), html.P("No. of variant observations")],
                                    className="mini_container",
                                ),
                            ],
                            id="info-container",
                            className="row container-display",
                        ),
                        html.Div(
                            [dcc.Graph(id="count_graph")],
                            id="countGraphContainer",
                            className="pretty_container",
                        ),
                    ],
                    id="right-column",
                    className="eight columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                html.Div(
                    [dcc.Graph(id="main_graph")],
                    className="pretty_container seven columns",
                ),
                html.Div(
                    [dcc.Graph(id="individual_graph")],
                    className="pretty_container five columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                html.Div(
                    [dcc.Graph(id="pie_graph")],
                    className="pretty_container seven columns",
                ),
                html.Div(
                    [dcc.Graph(id="aggregate_graph")],
                    className="pretty_container five columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Footer(
            [
                html.Div(
                    [
                        html.P(),
                        html.P("© 2021 TRON Mainz. All Rights Reserved"),
                        html.P()
                    ],
                    className="one-third column"
                ),
                html.Div(
                    [
                        html.P(),
                        html.P([html.A("DATA PROTECTION", href="https://tron-mainz.de/data-protection/"), " | ",
                                html.A("IMPRINT", href="https://tron-mainz.de/imprint/")]),
                        html.P()
                    ],
                    className="one-third column"
                )
            ]
        )
    ],
    id="mainContainer",
    style={"display": "flex", "flex-direction": "column"},
)


# Main graph -> individual graph
@app.callback(Output("individual_graph", "figure"), Input('individual_graph', 'clickData'),)
def stacked_area_samples_by_country(clickData):
    fig = px.area(data2_merged, x="date", y="cumsum", color="country",  # line_group="country",
                  category_orders={
                      "country": list(data2_merged.sort_values("cumsum", ascending=False).country.unique())[::-1]},
                  labels={"cumsum": "num. samples", "count": "increment"},
                  title="Accumulated samples per country",
                  hover_data=["count"],
                  color_discrete_sequence=random.shuffle(px.colors.qualitative.Dark24))
    fig.update_layout(
        legend={'traceorder': 'reversed'},
        xaxis={'title': None},
        yaxis={'dtick': 2000}
    )
    fig.add_trace(go.Scatter(
        x=["2021-01-17", "2021-01-17", "2021-01-17"],
        y=[16000, 9000, 2000],
        mode="text",
        text=["USA", "GBR", "IND"],
        textposition="top center",
        showlegend=False
    ))
    return fig


#@app.callback(
#    Output('example-graph', 'figure'),
#    Output('count_by_country', 'figure'),
#    Input('dropdown', 'value'))
def update_figure(selected_country):

    # filters data
    filtered_df = runs[runs.country.isin(selected_country)] \
        if selected_country is not None and len(selected_country) > 0 else runs

    fig = px.scatter(filtered_df, x="collection_date", y="country",
                     color="library_strategy",  hover_name="run_accession", opacity=0.8,
                     title="Samples through time") #, symbol="host_sex" size="read_count",)
                        # log_x=True, size_max=60)
    fig.update_layout(transition_duration=500)
    fig.update_layout(clickmode='event+select')
    fig.update_traces(marker_size=8)

    counts = filtered_df.groupby(["country", "library_strategy"]).size().reset_index(name="counts")
    fig2 = px.bar(counts, x="country", y="counts", color="library_strategy", title="Counts of samples")
    fig2.update_layout(xaxis_tickangle=-45)

    return fig, fig2


#@app.callback(
#    [Output("table", "data"), Output('table', 'columns')],
#    Input('example-graph', 'clickData'),
#    Input('example-graph', 'selectedData'),
#)
def display_click_data(clickData, selectedData):
    run_accession = None
    if clickData:
        points = clickData.get("points")
        if points and len(points) > 0:
            run_accession = [points[0].get("hovertext")]
    elif selectedData:
        points = selectedData.get("points")
        if points and len(points) > 0:
            run_accession = [p.get("hovertext") for p in points]

    data = runs
    if run_accession:
        data = runs[runs.run_accession.isin(run_accession)]
    return data.to_dict('records'), [{"name": i, "id": i} for i in runs.columns]



def main(debug=False):
    app.run_server(debug=debug)


if __name__ == '__main__':
    main(debug=True)
