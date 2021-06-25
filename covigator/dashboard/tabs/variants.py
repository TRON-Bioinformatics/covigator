from datetime import timedelta, datetime
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from dash.dependencies import Output, Input
from covigator.dashboard.figures import Figures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE, MONTH_PATTERN
from covigator.database.queries import Queries


def get_tab_variants(queries: Queries):

    return dcc.Tab(label="Variants",
                   style=TAB_STYLE,
                   selected_style=TAB_SELECTED_STYLE,
                   children=[
                       html.Div(
                           id='needleplot-body',
                           className="row container-display",
                           style={'overflow': 'scroll'}, # 'top': 0, 'bottom': 0, position: fixed
                           children=[
                               get_variants_tab_left_bar(queries=queries),
                               get_variants_tab_graphs()
                           ]),
                   ])


def get_variants_tab_graphs():
    return html.Div(children=[
        html.Br(),
        html.Div(id='top-occurring-variants', children=dash_table.DataTable(id="top-occurring-variants-table")),
        html.Br(),
        html.Div(id='needle-plot'),
        html.Br(),
        html.Div(id='cooccurrence-heatmap'),
        html.Br(),
        html.Div(id='variants-mds'),
        html.Br(),
    ], className="ten columns", style={'overflow': 'scroll', "height": "900px"}, )


def get_variants_tab_left_bar(queries: Queries):

    genes = queries.get_genes()
    months = queries.get_sample_months(MONTH_PATTERN)
    today = datetime.now()
    today_formatted = today.strftime(MONTH_PATTERN)
    oneyearago = today - timedelta(days=356)
    oneyearago_formatted = oneyearago.strftime(MONTH_PATTERN)

    return html.Div(children=[
        html.Br(),
        dcc.Markdown("Select a gene"),
        dcc.Dropdown(
            id='dropdown-gene',
            options=[{'label': c.name, 'value': c.name} for c in genes],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""**Top occurring variants**
        Number of top occurring variants"""),
        dcc.Slider(
            id='slider-top-variants',
            min=10,
            max=100,
            step=10,
            value=10,
            dots=True,
            marks={i: '{}'.format(i) for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]},
            tooltip=dict(always_visible=False, placement="right")
        ),
        html.Br(),
        dcc.Markdown("""Metric to measure abundance of variants per month"""),
        dcc.Dropdown(
            id='dropdown-top-variants-metric',
            options=[
                {'label': 'Count', 'value': 'count'},
                {'label': 'Frequency', 'value': 'frequency_by_month'}],
            value='count',
            clearable=False,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""Select a start and end date"""),
        html.Div(children=[
            html.Div(
                dcc.Dropdown(
                    id='dropdown-date-range-start',
                    options=[{'label': c, 'value': c} for c in months],
                    value=oneyearago_formatted,
                    multi=False,
                    clearable=False
                ), className="six column"),
            html.Div(id="dropdown-date-range-end-div",
                     children=dcc.Dropdown(
                         id='dropdown-date-range-end',
                         options=[{'label': c, 'value': c} for c in months],
                         value=today_formatted,
                         multi=False,
                         clearable=False
                     ),
                     className="six column"),
        ], className="row container-display"),
        html.Br(),
        dcc.Markdown("""**Genome view**
        
        Bin size"""),
        dcc.Slider(
            id='slider-bin-size',
            min=5,
            max=400,
            step=5,
            value=50,
            dots=False,
            marks={i: '{}'.format(i) for i in [10, 50, 100, 200, 300, 400]},
            tooltip=dict(always_visible=False, placement="right")
        ),
        html.Br(),
        dcc.Markdown("""**Co-occurrence matrix**
        
        Metric to assess paiwise co-occurrence"""),
        dcc.Dropdown(
            id='dropdown-heatmap-metric',
            options=[{'label': "Count", 'value': "count"},
                     {'label': "Frequency", 'value': "frequency"},
                     {'label': "Jaccard index", 'value': "jaccard"},
                     {'label': "Cohen's kappa coefficient", 'value': "kappa"},
                     ],
            value="jaccard",
            clearable=False,
            multi=False
        ),
        dcc.Markdown("""Minimum number of pairwise co-occurrences"""),
        dcc.Slider(
            id='slider-min-cooccurrences',
            min=1,
            max=100,
            step=5,
            value=10,
            dots=True,
            marks={i: '{}'.format(i) for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]},
            tooltip=dict(always_visible=False, placement="right")
        ),
        html.Br(),
        dcc.Markdown("""**Variants clustering**
        
        The number of samples (or total weight) in a neighborhood for a point to be 
        considered as a core point. This includes the point itself."""),
        dcc.Slider(
            id='slider-min-samples',
            min=2,
            max=10,
            step=1,
            value=5,
            dots=True,
            marks={i: str(i) for i in range(2, 10)},
            tooltip=dict(always_visible=False, placement="right")
        ),
    ], className="two columns")


def set_callbacks_variants_tab(app, figures: Figures, queries: Queries):
    @app.callback(
        Output('top-occurring-variants', 'children'),
        Input('slider-top-variants', 'value'),
        Input('dropdown-gene', 'value'),
        Input('dropdown-date-range-start', 'value'),
        Input('dropdown-date-range-end', 'value'),
        Input('dropdown-top-variants-metric', 'value')
    )
    def update_top_occurring_variants(top_variants, gene_name, date_range_start, date_range_end, metric):
        return html.Div(children=figures.get_top_occurring_variants_plot(
            top=top_variants, gene_name=gene_name, date_range_start=date_range_start,
            date_range_end=date_range_end, metric=metric))

    @app.callback(
        Output('needle-plot', 'children'),
        Input('dropdown-gene', 'value'),
        Input('top-occurring-variants-table', "derived_virtual_data"),
        Input('top-occurring-variants-table', "derived_virtual_selected_rows"),
        Input('slider-bin-size', 'value'),
    )
    def update_needle_plot(gene_name, rows, selected_rows_indices, bin_size):
        if gene_name is not None:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(
                children=figures.get_variants_plot(
                    gene_name=gene_name,
                    selected_variants=selected_rows,
                    bin_size=bin_size))
        else:
            plot = html.Div(
                children=figures.get_variants_abundance_plot(
                    bin_size=bin_size))
        return plot

    @app.callback(
        Output('dropdown-date-range-end-div', 'children'),
        Input('dropdown-date-range-start', 'value'))
    def update_dropdown_end_date(start_date):
        today = datetime.now()
        today_formatted = today.strftime(MONTH_PATTERN)
        months = [m for m in queries.get_sample_months(MONTH_PATTERN) if m >= start_date]
        return dcc.Dropdown(
            id='dropdown-date-range-end',
            options=[{'label': c, 'value': c} for c in months],
            value=today_formatted,
            multi=False,
            clearable=False
        )

    @app.callback(
        Output('cooccurrence-heatmap', 'children'),
        Input('dropdown-gene', 'value'),
        Input('top-occurring-variants-table', "derived_virtual_data"),
        Input('top-occurring-variants-table', "derived_virtual_selected_rows"),
        Input('dropdown-heatmap-metric', 'value'),
        Input('slider-min-cooccurrences', 'value'),
    )
    def update_cooccurrence_heatmap(gene_name, rows, selected_rows_indices, metric, min_occurrences):
        selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
        plot = html.Div(children=figures.get_cooccurrence_heatmap(
            gene_name=gene_name,
            selected_variants=selected_rows,
            metric=metric,
            min_occurrences=min_occurrences))
        return plot

    @app.callback(
        Output('variants-mds', 'children'),
        Input('dropdown-gene', 'value'),
        Input('top-occurring-variants-table', "derived_virtual_data"),
        Input('top-occurring-variants-table', "derived_virtual_selected_rows"),
        Input('slider-min-cooccurrences', 'value'),
        Input('slider-min-samples', 'value')
    )
    def update_variants_mds(gene_name, rows, selected_rows_indices, min_cooccurrence, min_samples):
        selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
        plot = html.Div(children=figures.get_variants_clustering(
            gene_name=gene_name,
            selected_variants=selected_rows,
            min_cooccurrence=min_cooccurrence,
            min_samples=min_samples))
        return plot