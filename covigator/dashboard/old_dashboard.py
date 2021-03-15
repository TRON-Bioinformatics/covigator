# Run this app with `python old_dashboard.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import plotly.express as px
import pandas as pd
from dash.dependencies import Output, Input
from covigator.model import SampleEna
from covigator.database import Database

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
session = Database().get_database_session()

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

# read all samples data
runs = pd.read_sql(session.query(SampleEna).statement, session.bind)

# sets up page layout

app.layout = html.Div(children=[
    html.H1(children='Covigator'),
    dcc.Markdown('''
    ### Covigator dash prototype
    This is just to show off what could be done with dash **without needing to setup an API or implement a
    dedicated frontend** https://dash.plotly.com. Dash is a dashboard based on plotly interactive visualizations.

    This text is written in markup language by the way http://commonmark.org/.
    
    This prototype only uses the table of runs, loads it from the database in memory and visualizations are updated with callback functions.
    For bigger datasets the callback functions could be triggering queries to the database.
    '''),
    html.Br(),
      html.Div(className='row',  # Define the row element
               children=[
                   html.Div(className='three columns div-user-controls', children=[
                        html.Div([
                             dcc.Markdown("""
                             #### Total samples: {}
                             #### Number of countries: {}
                             #### First sample: {}
                             #### Last sample: {}
                             """.format(runs.shape[0], len(set(runs.country)),
                                        min(runs.collection_date[runs.collection_date != ""].dropna()),
                                        max(runs.collection_date[runs.collection_date != ""].dropna())))
                        ]),
                        html.Br(),
                        html.Div([
                             dcc.Markdown("""
                             **Select one or more countries**
                             """),
                            dcc.Dropdown(
                                id='dropdown',
                                options=[{'label': c, 'value': c} for c in
                                         sorted(set(list(runs.country)))],
                                value=None,
                                multi=True
                            ),
                        ]),
                   ]),
                   html.Div(className='nine columns div-for-charts bg-grey', children=[
                       dcc.Graph(id='example-graph'),
                       dcc.Graph(id='count_by_country'),
                       ])
               ]),
      html.Div(className="row", children=[
          dash_table.DataTable(
              id='table',
              data=[],
              page_size=10,
            style_cell_conditional=[
                {
                    'if': {'column_id': c},
                    'textAlign': 'left'
                } for c in ['Date', 'Region']
            ],
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
            style_header={
                'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold'
            },
            sort_action='native',
            filter_action='native',
            hidden_columns=["scientific_name", "library_name", "nominal_lenght", "library_source", "library_selection",
                            "sample_collection", "sequencing_method", "fastq_ftp", "fastq_md5", "host_tax_id",
                            "sample_accession", "study_accession", "experiment_accession"]
              )
      ])

])


@app.callback(
    Output('example-graph', 'figure'),
    Output('count_by_country', 'figure'),
    Input('dropdown', 'value'))
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


@app.callback(
    [Output("table", "data"), Output('table', 'columns')],
    Input('example-graph', 'clickData'),
    Input('example-graph', 'selectedData'),
)
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
