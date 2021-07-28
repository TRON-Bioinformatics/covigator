import functools
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
from dash.dependencies import Output, Input
from sqlalchemy.orm import Session

from covigator.dashboard.figures.subclonal_variants import SubclonalVariantsFigures
from covigator.database.queries import Queries

TOP_COOCCURRING_CLONAL_VARIANTS = 'top-cooccurring-clonal-variants-table'

ID_HIST_LIBRARY_STRATEGY = "id-hist-library-strategy"
ID_HIST_COUNTRIES = "id-hist-countries"

ID_DROPDOWN_ORDER_BY = "id-combobox-order-by"

ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS = 'dropdown-gene-subclonal-variants'
ID_SLIDER_SUBCLONAL_VARIANTS_VAF = 'slider-subclonal-variants-vaf'
ID_SLIDER_TOP_SUBCLONAL_VARIANTS = 'slider-subclonal-variants-top'

ID_TOP_OCCURRING_SUBCLONAL_VARIANTS = 'top-occurring-subclonal-variants'
ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE = 'top-occurring-subclonal-variants-table'


@functools.lru_cache()
def get_tab_subclonal_variants(queries: Queries):

    return dbc.Card(
        dbc.CardBody(
            children=[
                get_subclonal_variants_tab_left_bar(queries=queries),
                get_subclonal_variants_tab_graphs()
            ])
    )


def get_subclonal_variants_tab_graphs():
    return html.Div(children=[
        html.Br(),
        html.Div(id=ID_TOP_OCCURRING_SUBCLONAL_VARIANTS,
                 children=dash_table.DataTable(id=ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE)),
        html.Br(),
        html.Div(id=ID_HIST_LIBRARY_STRATEGY,
                 className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
        html.Div(id=ID_HIST_COUNTRIES,
                 className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
        html.Br(),
        html.Div(id=TOP_COOCCURRING_CLONAL_VARIANTS,
                 className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
    ], className="ten columns", style={'overflow': 'scroll', "height": "900px"}, )


def get_subclonal_variants_tab_left_bar(queries: Queries):

    genes = queries.get_genes()

    return html.Div(children=[
        html.Br(),
        dcc.Markdown("Select a gene"),
        dcc.Dropdown(
            id=ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS,
            options=[{'label': c.name, 'value': c.name} for c in genes],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""Minimum VAF subclonal variants"""),
        dcc.Slider(
            id=ID_SLIDER_SUBCLONAL_VARIANTS_VAF,
            min=0.1,
            max=0.8,
            step=0.05,
            value=0.3,
            marks={i: '{}'.format(i) for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]},
            tooltip=dict(always_visible=False, placement="right")
        ),
        html.Br(),
        dcc.Markdown("""Number of top occurring subclonal variants"""),
        dcc.Slider(
            id=ID_SLIDER_TOP_SUBCLONAL_VARIANTS,
            min=10,
            max=200,
            step=10,
            value=10,
            marks={i: '{}'.format(i) for i in [10, 50, 100, 150, 200]},
            tooltip=dict(always_visible=False, placement="right")
        ),
        html.Br(),
        dcc.Markdown("Order results by"),
        dcc.Dropdown(
            id=ID_DROPDOWN_ORDER_BY,
            options=[{'label': 'Score (ConsHMM + count)', 'value': 'score'},
                     {'label': 'Count observations', 'value': 'count'},
                     {'label': 'ConsHMM', 'value': 'conservation'},
                     {'label': 'VAF', 'value': 'vaf'}],
            value='score',
            multi=False
        ),
    ], className="two columns")


def set_callbacks_subclonal_variants_tab(app, session: Session):

    queries = Queries(session=session)
    figures = SubclonalVariantsFigures(queries=queries)

    @app.callback(
        Output(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS, 'children'),
        Input(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
        Input(ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS, 'value'),
        Input(ID_SLIDER_TOP_SUBCLONAL_VARIANTS, 'value'),
        Input(ID_DROPDOWN_ORDER_BY, 'value')
    )
    def update_top_occurring_variants(min_vaf, gene_name, top, order_by):
        return html.Div(children=figures.get_top_occurring_subclonal_variants_plot(
            min_vaf=min_vaf, gene_name=gene_name, top=top, order_by=order_by))

    @app.callback(
        Output(ID_HIST_LIBRARY_STRATEGY, 'children'),
        Input(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
    )
    def update_hist_library_strategy(min_vaf, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_hist_library_strategy(variant_id=variant_id, min_vaf=min_vaf),
            )
        return plot

    @app.callback(
        Output(ID_HIST_COUNTRIES, 'children'),
        Input(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
    )
    def update_hist_countries(min_vaf, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_hist_countries(variant_id=variant_id, min_vaf=min_vaf),
            )
        return plot

    @app.callback(
        Output(TOP_COOCCURRING_CLONAL_VARIANTS, 'children'),
        Input(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
        Input(ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS, 'value'),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
        Input(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
    )
    def update_cooccurring_clonal_variants(min_vaf, gene_name, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_cooccurring_clonal_variants(variant_id=variant_id, min_vaf=min_vaf,
                                                                 gene_name=gene_name),
            )
        return plot
