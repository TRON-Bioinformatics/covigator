from datetime import datetime

from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash import dash_table
from covigator.dashboard.tabs import get_mini_container, print_number
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session
import pandas as pd

from covigator.dashboard.figures.intrahost_mutations import IntrahostMutationsFigures
from covigator.database.queries import Queries

TOP_COOCCURRING_CLONAL_VARIANTS = 'top-cooccurring-clonal-variants-table'

ID_HIST_LIBRARY_STRATEGY = "id-hist-library-strategy"
ID_HIST_COUNTRIES = "id-hist-countries"

ID_DROPDOWN_ORDER_BY = "id-combobox-order-by"

ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS = 'dropdown-gene-subclonal-variants'
ID_DROPDOWN_DOMAIN_SUBCLONAL_VARIANTS = 'dropdown-domain-subclonal-variants'
ID_SLIDER_SUBCLONAL_VARIANTS_VAF = 'slider-subclonal-variants-vaf'
ID_SLIDER_TOP_SUBCLONAL_VARIANTS = 'slider-subclonal-variants-top'

ID_TOP_OCCURRING_SUBCLONAL_VARIANTS = 'top-occurring-subclonal-variants'
ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE = 'top-occurring-subclonal-variants-table'
ID_APPLY_BUTTOM = 'im-apply-buttom'


def get_tab_subclonal_variants(queries: Queries):

    return dbc.CardBody(
            children=[
                get_subclonal_variants_tab_left_bar(queries=queries),
                html.Div(
                    className="one columns",
                    children=[html.Br()]),
                get_subclonal_variants_tab_graphs(queries)
            ])


def get_subclonal_variants_tab_graphs(queries):

    return html.Div(
        children=[
            html.Div(id=ID_TOP_OCCURRING_SUBCLONAL_VARIANTS,
                     children=dash_table.DataTable(id=ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE)),
            html.Hr(),
            html.Br(),
            html.Div(id=ID_HIST_LIBRARY_STRATEGY,
                     className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
            html.Div(id=ID_HIST_COUNTRIES,
                     className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
            html.Hr(),
            html.Br(),
            html.Div(id=TOP_COOCCURRING_CLONAL_VARIANTS,
                     className="five columns", style={"margin-left": 0, "margin-right": "1%", "width": "48%"}),
        ],
        className="nine columns",
    )


def get_subclonal_variants_tab_left_bar(queries: Queries):

    genes = queries.get_genes()
    count_subclonal_variant_observations = queries.count_subclonal_variant_observations()
    count_unique_subclonal_variant = queries.count_unique_subclonal_variant()
    count_unique_only_subclonal_variant = queries.count_unique_only_subclonal_variant()

    return html.Div(children=[
        dcc.Markdown("""
                    Mutations with a VAF within the interval [0.05, 0.8) are considered intrahost variants. 
                    Intrahost mutations can only be detected when the raw reads are available.
                    The dataset of intrahost mutations is enriched for false positive calls due to the lower 
                    Variant Allele Frequency (VAF). (Lythgoe, 2021) and (Valesano, 2021) reported that SARS-CoV-2 
                    intrahost variant calls with a VAF below 0.02 and 0.03 respectively had not enough quality; 
                    although the variant calling methods differ between them and with CoVigator.

                    Here we provide a tool to explore most frequent intrahost mutations that have not been observed 
                    before as clonal variants.
                     """, style={"font-size": 16}),
        html.Br(),
        html.Div(
            html.Span(
                children=[
                    get_mini_container(
                        title="Unique mutations",
                        value=print_number(count_unique_subclonal_variant)
                    ),
                    html.Br(),
                    html.Br(),
                    get_mini_container(
                        title="Only intrahost",
                        value=print_number(count_unique_only_subclonal_variant)
                    ),
                    html.Br(),
                    html.Br(),
                    get_mini_container(
                        title="Mutation calls",
                        value=print_number(count_subclonal_variant_observations)
                    ),
                    html.Br(),
                ]
            )
        ),
        html.Br(),
        dcc.Markdown("Select a gene"),
        dcc.Dropdown(
            id=ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS,
            options=[{'label': c.name, 'value': c.name} for c in genes],
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""Select a protein domain"""),
        dcc.Dropdown(
            id=ID_DROPDOWN_DOMAIN_SUBCLONAL_VARIANTS,
            value=None,
            multi=False
        ),
        html.Br(),
        dcc.Markdown("""Minimum VAF intrahost variants"""),
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
        dcc.Markdown("""Number of top occurring intrahost variants"""),
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
        html.Br(),
        html.Button('Apply', id=ID_APPLY_BUTTOM),
    ], className="two columns")


def set_callbacks_subclonal_variants_tab(app, session: Session):

    queries = Queries(session=session)
    figures = IntrahostMutationsFigures(queries=queries)
    domains_by_gene = {g.name: queries.get_domains_by_gene(g.name) for g in queries.get_genes()}
    all_domains = queries.get_domains()

    @app.callback(
        Output(ID_DROPDOWN_DOMAIN_SUBCLONAL_VARIANTS, 'options'),
        Input(ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS, 'value'))
    def set_domains(selected_gene):
        domains = domains_by_gene.get(selected_gene) if selected_gene else all_domains
        domain_labels = sorted({("{gene}: {domain}".format(domain=d.name, gene=d.gene_name), d.name) for d in domains})
        return [{'label': label, 'value': value} for label, value in domain_labels]

    @app.callback(
        Output(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
            State(ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS, 'value'),
            State(ID_DROPDOWN_DOMAIN_SUBCLONAL_VARIANTS, 'value'),
            State(ID_SLIDER_TOP_SUBCLONAL_VARIANTS, 'value'),
            State(ID_DROPDOWN_ORDER_BY, 'value')
        ]
    )
    def update_top_occurring_variants(_, min_vaf, gene_name, domain, top, order_by):
        return html.Div(children=figures.get_top_occurring_subclonal_variants_plot(
            min_vaf=min_vaf, gene_name=gene_name, domain=domain, top=top, order_by=order_by))

    @app.callback(
        Output(ID_HIST_LIBRARY_STRATEGY, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
        ]
    )
    def update_hist_library_strategy(_, min_vaf, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_hist_library_strategy(variant_id=variant_id, min_vaf=min_vaf),
            )
        return plot

    @app.callback(
        Output(ID_HIST_COUNTRIES, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
        ]
    )
    def update_hist_countries(_, min_vaf, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_hist_countries(variant_id=variant_id, min_vaf=min_vaf),
            )
        return plot

    @app.callback(
        Output(TOP_COOCCURRING_CLONAL_VARIANTS, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_SLIDER_SUBCLONAL_VARIANTS_VAF, 'value'),
            State(ID_DROPDOWN_GENE_SUBCLONAL_VARIANTS, 'value'),
            State(ID_DROPDOWN_DOMAIN_SUBCLONAL_VARIANTS, 'value'),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_data"),
            State(ID_TOP_OCCURRING_SUBCLONAL_VARIANTS_TABLE, "derived_virtual_selected_rows")
        ]
    )
    def update_cooccurring_clonal_variants(_, min_vaf, gene_name, domain, data, selected_rows):
        plot = None
        if selected_rows:
            variant_id = data[selected_rows[0]].get("variant_id")
            plot = html.Div(
                children=figures.get_cooccurring_clonal_variants(
                    variant_id=variant_id, min_vaf=min_vaf, gene_name=gene_name, domain=domain),
            )
        return plot

    @app.callback(
        Output("download-dataframe-csv3", "data"),
        inputs=[
            Input("btn_csv3", "n_clicks"),
            Input("memory3", "data")
        ],
        prevent_initial_call=True,
    )
    def download_top_intrahost_mutations(n_clicks, df):
        return dcc.send_data_frame(
            pd.DataFrame.from_dict(df).to_csv,
            "covigator_top_intrahost_mutations_{}.csv".format(datetime.now().strftime("%Y%m%d%H%M%S")),
            index=False
        )
