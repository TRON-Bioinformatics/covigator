import functools
from datetime import timedelta, datetime
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session
from covigator.dashboard.figures.recurrent_mutations import RecurrentMutationsFigures
from covigator.dashboard.tabs import MONTH_PATTERN
from covigator.database.model import DataSource
from covigator.database.queries import Queries

ID_DROPDOWN_DATE_RANGE_END_DIV = 'dropdown-date-range-end-div'
ID_DROPDOWN_GENE = 'dropdown-gene'
ID_DROPDOWN_DOMAIN = 'dropdown-domain'
ID_SLIDER_MIN_SAMPLES = 'slider-min-samples'
ID_SLIDER_MIN_COOCCURRENCES = 'slider-min-cooccurrences'
ID_DROPDOWN_SIMILARITY_METRIC = 'dropdown-heatmap-metric'
ID_SLIDER_BIN_SIZE = 'slider-bin-size'
BIN_SIZE_VALUES = [10, 50, 100, 200, 300, 400]
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
ID_APPLY_BUTTOM = 'rm-apply-buttom'



@functools.lru_cache()
def get_tab_variants(queries: Queries, data_source: DataSource):

    return dbc.CardBody(
            children=[
                get_variants_tab_left_bar(queries=queries, data_source=data_source),
                html.Div(
                    className="one columns",
                    children=[html.Br()]),
                get_variants_tab_graphs()
            ])


def get_variants_tab_graphs():
    return html.Div(
        children=[
            html.Br(),
            html.Div(id=ID_TOP_OCCURRING_VARIANTS, children=dash_table.DataTable(id=ID_TOP_OCCURRING_VARIANTS_TABLE)),
            html.Br(),
            html.Div(id=ID_NEEDLE_PLOT),
            html.Br(),
            html.Div(id=ID_COOCCURRENCE_HEATMAP),
            html.Br(),
            html.Div(id=ID_VARIANTS_MDS),
            html.Br(),
        ],
        className="nine columns",
    )


def get_variants_tab_left_bar(queries: Queries, data_source: DataSource):

    # removes repeated gene names (ie: ORF1ab)
    genes = sorted({c.name for c in queries.get_genes()})
    months = queries.get_sample_months(MONTH_PATTERN)
    today = datetime.now()
    today_formatted = today.strftime(MONTH_PATTERN)
    oneyearago = today - timedelta(days=356)
    oneyearago_formatted = oneyearago.strftime(MONTH_PATTERN)

    return html.Div(children=[
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
        dcc.Markdown("Select a gene"),
        dcc.Dropdown(
            id=ID_DROPDOWN_GENE,
            options=[{'label': g, 'value': g} for g in genes],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""Select a protein domain"""),
        dcc.Dropdown(
            id=ID_DROPDOWN_DOMAIN,
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""**Top occurring mutations table**

Number of top occurring mutations"""),
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
        dcc.Markdown("""Metric to measure abundance of mutations per month"""),
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
            dcc.Dropdown(
                id=ID_DROPDOWN_DATE_RANGE_START,
                options=[{'label': c, 'value': c} for c in months],
                value=oneyearago_formatted,
                multi=False,
                clearable=False
            ),
            html.Div(
                id=ID_DROPDOWN_DATE_RANGE_END_DIV,
                children=dcc.Dropdown(
                    id=ID_DROPDOWN_DATE_RANGE_END,
                    options=[{'label': c, 'value': c} for c in months],
                    value=today_formatted,
                    multi=False,
                    clearable=False
                ))]),
        html.Br(),
        dcc.Markdown("""**Genome view**
        
Bin size"""),
        dcc.Dropdown(
            id=ID_SLIDER_BIN_SIZE,
            options=[{'label': "{} bp".format(c), 'value': str(c)} for c in BIN_SIZE_VALUES],
            value=50,
            multi=False,
            clearable=False
        ),
        html.Br(),
        dcc.Markdown("""**Co-occurrence matrix**
        
Paiwise co-occurrence metric"""),
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
        dcc.Markdown("""Minimum pairwise co-occurrences"""),
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
        dcc.Markdown("""**Mutations clustering**
        
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
        html.Br(),
        html.Button('Apply', id=ID_APPLY_BUTTOM),
    ], className="two columns")


def set_callbacks_variants_tab(app, session: Session):

    queries = Queries(session=session)
    figures = RecurrentMutationsFigures(queries=queries)
    domains_by_gene = {g.name: queries.get_domains_by_gene(g.name) for g in queries.get_genes()}
    all_domains = queries.get_domains()
    months_from_db = queries.get_sample_months(MONTH_PATTERN)

    @app.callback(
        Output(ID_DROPDOWN_DOMAIN, 'options'),
        Input(ID_DROPDOWN_GENE, 'value'))
    def set_domains(selected_gene):
        domains = domains_by_gene.get(selected_gene) if selected_gene else all_domains
        domain_labels = sorted({("{gene}: {domain}".format(domain=d.name, gene=d.gene_name), d.name) for d in domains})
        return [{'label': label, 'value': value} for label, value in domain_labels]

    @app.callback(
        Output(ID_DROPDOWN_DATE_RANGE_END_DIV, 'children'),
        Input(ID_DROPDOWN_DATE_RANGE_START, 'value'),
    )
    def update_dropdown_end_date(start_date):
        today = datetime.now()
        today_formatted = today.strftime(MONTH_PATTERN)
        months = [m for m in months_from_db if m >= start_date]
        return dcc.Dropdown(
            id=ID_DROPDOWN_DATE_RANGE_END,
            options=[{'label': c, 'value': c} for c in months],
            value=today_formatted,
            multi=False,
            clearable=False
        )

    @app.callback(
        Output(ID_TOP_OCCURRING_VARIANTS, 'children'),
        [Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_SLIDER_TOP_VARIANTS, 'value'),
            State(ID_DROPDOWN_GENE, 'value'),
            State(ID_DROPDOWN_DOMAIN, 'value'),
            State(ID_DROPDOWN_DATE_RANGE_START, 'value'),
            State(ID_DROPDOWN_DATE_RANGE_END, 'value'),
            State(ID_TOP_VARIANTS_METRIC, 'value'),
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
        ]
    )
    def update_top_occurring_variants(_, top_variants, gene_name, domain, date_range_start, date_range_end, metric, source):
        return html.Div(children=figures.get_top_occurring_variants_plot(
            top=top_variants, gene_name=gene_name, domain=domain, date_range_start=date_range_start,
            date_range_end=date_range_end, metric=metric, source=source))

    @app.callback(
        Output(ID_NEEDLE_PLOT, 'children'),
        [Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_GENE, 'value'),
            State(ID_DROPDOWN_DOMAIN, 'value'),
            State(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_data"),
            State(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_selected_rows"),
            State(ID_SLIDER_BIN_SIZE, 'value'),
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
        ]
    )
    def update_needle_plot(_, gene_name, domain, rows, selected_rows_indices, bin_size, source):
        if gene_name is not None or domain is not None:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(
                children=figures.get_variants_plot(
                    gene_name=gene_name,
                    domain_name=domain,
                    selected_variants=selected_rows,
                    bin_size=int(bin_size),
                    source=source
                ))
        else:
            plot = html.Div(
                children=figures.get_variants_abundance_plot(bin_size=int(bin_size), source=source))
        return plot

    @app.callback(
        Output(ID_COOCCURRENCE_HEATMAP, 'children'),
        [Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_GENE, 'value'),
            State(ID_DROPDOWN_DOMAIN, 'value'),
            State(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_data"),
            State(ID_TOP_OCCURRING_VARIANTS_TABLE, "derived_virtual_selected_rows"),
            State(ID_DROPDOWN_SIMILARITY_METRIC, 'value'),
            State(ID_SLIDER_MIN_COOCCURRENCES, 'value'),
            State(ID_SLIDER_MIN_SAMPLES, 'value'),
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
        ]
    )
    def update_cooccurrence_heatmap(
            _, gene_name, domain, rows, selected_rows_indices, metric, min_cooccurrences, min_samples, source):
        if source != DataSource.ENA.name:
            plot = html.Div(
                children=[dcc.Markdown(
                    """**The co-occurrence analysis is currently only available for the ENA dataset**""")])
        elif gene_name is None:
            plot = html.Div(
                children=[dcc.Markdown(
                    """**Please, select a gene or domain to explore the co-occurrence analysis**""")])
        else:
            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            sparse_matrix = queries.get_sparse_cooccurrence_matrix(
                gene_name=gene_name, domain=domain, min_cooccurrence=min_cooccurrences)
            plot = html.Div(children=[
                figures.get_cooccurrence_heatmap(
                    sparse_matrix=sparse_matrix,
                    selected_variants=selected_rows,
                    metric=metric,
                    min_cooccurrences=min_cooccurrences),
                figures.get_variants_clustering(
                    sparse_matrix=sparse_matrix, min_cooccurrence=min_cooccurrences, min_samples=min_samples)
            ]
            )
        return plot
