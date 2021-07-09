from datetime import timedelta, datetime
import dash_core_components as dcc
import dash_html_components as html
import dash_table
from dash.dependencies import Output, Input
from sqlalchemy.orm import Session
from covigator.dashboard.figures.variants import VariantsFigures
from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE, MONTH_PATTERN
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_DROPDOWN_DATE_RANGE_END_DIV = 'dropdown-date-range-end-div'
ID_DROPDOWN_GENE = 'dropdown-gene'
ID_SLIDER_MIN_SAMPLES = 'slider-min-samples'
ID_SLIDER_MIN_COOCCURRENCES = 'slider-min-cooccurrences'
ID_DROPDOWN_SIMILARITY_METRIC = 'dropdown-heatmap-metric'
ID_SLIDER_BIN_SIZE = 'slider-bin-size'
ID_DROPDOWN_DATE_RANGE_END = 'dropdown-date-range-end'
ID_DROPDOWN_DATE_RANGE_START = 'dropdown-date-range-start'
ID_TOP_VARIANTS_METRIC = 'dropdown-top-variants-metric'
ID_SLIDER_TOP_VARIANTS = 'slider-top-variants'
ID_VARIANTS_MDS = 'variants-mds'
ID_COOCCURRENCE_HEATMAP = 'cooccurrence-heatmap'
ID_NEEDLE_PLOT = 'needle-plot'
ID_TOP_OCCURRING_VARIANTS = 'top-occurring-variants'
ID_TOP_OCCURRING_VARIANTS_TABLE = 'top-occurring-variants-table'
ID_DROPDOWN_DATA_SOURCE = "dropdown-data-source-variants-tab"


def get_tab_variants(queries: Queries):

    return dcc.Tab(label="Recurrent mutations",
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
        html.Div(id=ID_TOP_OCCURRING_VARIANTS, children=dash_table.DataTable(id=ID_TOP_OCCURRING_VARIANTS_TABLE)),
        html.Br(),
        html.Div(id=ID_NEEDLE_PLOT),
        html.Br(),
        html.Div(id=ID_COOCCURRENCE_HEATMAP),
        html.Br(),
        html.Div(id=ID_VARIANTS_MDS),
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
        dcc.Markdown("""Select data source"""),
        dcc.Dropdown(
            id=ID_DROPDOWN_DATA_SOURCE,
            options=[{'label': DataSource.ENA.name, 'value': DataSource.ENA.name},
                     {'label': DataSource.GISAID.name, 'value': DataSource.GISAID.name}],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("Select a gene"),
        dcc.Dropdown(
            id=ID_DROPDOWN_GENE,
            options=[{'label': c.name, 'value': c.name} for c in genes],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""**Top occurring variants**

Number of top occurring variants"""),
        dcc.Slider(
            id=ID_SLIDER_TOP_VARIANTS,
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
            id=ID_TOP_VARIANTS_METRIC,
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
                    id=ID_DROPDOWN_DATE_RANGE_START,
                    options=[{'label': c, 'value': c} for c in months],
                    value=oneyearago_formatted,
                    multi=False,
                    clearable=False
                ), className="six column"),
            html.Div(
                id=ID_DROPDOWN_DATE_RANGE_END_DIV,
                children=dcc.Dropdown(
                    id=ID_DROPDOWN_DATE_RANGE_END,
                    options=[{'label': c, 'value': c} for c in months],
                    value=today_formatted,
                    multi=False,
                    clearable=False
                ), className="six column"),
        ], className="row container-display"),
        html.Br(),
        dcc.Markdown("""**Genome view**
        
Bin size"""),
        dcc.Slider(
            id=ID_SLIDER_BIN_SIZE,
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
            id=ID_DROPDOWN_SIMILARITY_METRIC,
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
            id=ID_SLIDER_MIN_COOCCURRENCES,
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
            id=ID_SLIDER_MIN_SAMPLES,
            min=2,
            max=10,
            step=1,
            value=5,
            dots=True,
            marks={i: str(i) for i in range(2, 10)},
            tooltip=dict(always_visible=False, placement="right")
        ),
    ], className="two columns")


def set_callbacks_variants_tab(app, session: Session):

    queries = Queries(session=session)
    figures = VariantsFigures(queries=queries)

    @app.callback(
        Output(ID_TOP_OCCURRING_VARIANTS, 'children'),
        Input(ID_SLIDER_TOP_VARIANTS, 'value'),
        Input(ID_DROPDOWN_GENE, 'value'),
        Input(ID_DROPDOWN_DATE_RANGE_START, 'value'),
        Input(ID_DROPDOWN_DATE_RANGE_END, 'value'),
        Input(ID_TOP_VARIANTS_METRIC, 'value'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value')
    )
    def update_top_occurring_variants(top_variants, gene_name, date_range_start, date_range_end, metric, source):
        return html.Div(children=figures.get_top_occurring_variants_plot(
            top=top_variants, gene_name=gene_name, date_range_start=date_range_start,
            date_range_end=date_range_end, metric=metric, source=source))

    @app.callback(
        Output(ID_NEEDLE_PLOT, 'children'),
        Input(ID_DROPDOWN_GENE, 'value'),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_selected_rows"),
        Input(ID_SLIDER_BIN_SIZE, 'value'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value')
    )
    def update_needle_plot(gene_name, rows, selected_rows_indices, bin_size, source):
        if gene_name is not None:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(
                children=figures.get_variants_plot(
                    gene_name=gene_name,
                    selected_variants=selected_rows,
                    bin_size=bin_size,
                    source=source
                ))
        else:
            plot = html.Div(
                children=figures.get_variants_abundance_plot(bin_size=bin_size, source=source))
        return plot

    @app.callback(
        Output(ID_DROPDOWN_DATE_RANGE_END_DIV, 'children'),
        Input(ID_DROPDOWN_DATE_RANGE_START, 'value'))
    def update_dropdown_end_date(start_date):
        today = datetime.now()
        today_formatted = today.strftime(MONTH_PATTERN)
        months = [m for m in queries.get_sample_months(MONTH_PATTERN) if m >= start_date]
        return dcc.Dropdown(
            id=ID_DROPDOWN_DATE_RANGE_END,
            options=[{'label': c, 'value': c} for c in months],
            value=today_formatted,
            multi=False,
            clearable=False
        )

    @app.callback(
        Output(ID_COOCCURRENCE_HEATMAP, 'children'),
        Input(ID_DROPDOWN_GENE, 'value'),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_selected_rows"),
        Input(ID_DROPDOWN_SIMILARITY_METRIC, 'value'),
        Input(ID_SLIDER_MIN_COOCCURRENCES, 'value'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value')
    )
    def update_cooccurrence_heatmap(gene_name, rows, selected_rows_indices, metric, min_occurrences, source):
        if source != DataSource.ENA.name:
            plot = html.Div(
                children=[dcc.Markdown(
                    """**The cooccurrence heatmap is currently only available for the ENA dataset**""")])
        else:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(children=figures.get_cooccurrence_heatmap(
                gene_name=gene_name,
                selected_variants=selected_rows,
                metric=metric,
                min_occurrences=min_occurrences))
        return plot

    @app.callback(
        Output(ID_VARIANTS_MDS, 'children'),
        Input(ID_DROPDOWN_GENE, 'value'),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_selected_rows"),
        Input(ID_SLIDER_MIN_COOCCURRENCES, 'value'),
        Input(ID_SLIDER_MIN_SAMPLES, 'value'),
        Input(ID_DROPDOWN_DATA_SOURCE, 'value')
    )
    def update_variants_mds(gene_name, rows, selected_rows_indices, min_cooccurrence, min_samples, source):
        if source != DataSource.ENA.name:
            plot = html.Div(children=
                            [dcc.Markdown(
                                """**The variants clustering is currently only available for the ENA dataset**""")])
        else:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(children=figures.get_variants_clustering(
                gene_name=gene_name,
                selected_variants=selected_rows,
                min_cooccurrence=min_cooccurrence,
                min_samples=min_samples))
        return plot