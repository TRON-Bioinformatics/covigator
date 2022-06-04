from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Output, Input, State
from sqlalchemy.orm import Session

from covigator.dashboard.figures.mutation_stats import MutationStatsFigures
from covigator.dashboard.figures.samples import SampleFigures
from covigator.database.model import DataSource, VariantType
from covigator.database.queries import Queries


ID_VARIANTS_PER_SAMPLE_GRAPH = 'ms-variants-per-sample-graph'
ID_INDEL_LENGTH_GRAPH = 'ms-indel-lengths-graph'
ID_ANNOTATIONS_GRAPH = 'ms-id-annotations-graph'
ID_SUBSTITUTIONS_GRAPH = 'ms-substitutions-graph'
ID_SLIDER_MIN_SAMPLES = 'ms-slider-min-samples-per-country'
ID_DROPDOWN_DATA_SOURCE = "ms-dropdown-data-source"
ID_DROPDOWN_COUNTRY = 'ms-dropdown-country'
ID_DROPDOWN_GENE = 'ms-dropdown-gene-overall-mutations'
ID_DROPDOWN_VARIANT_TYPE = 'ms-dropdown-variant-type'
ID_ACCUMULATED_SAMPLES_GRAPH = 'ms-accumulated-samples-per-country'
ID_DN_DS_GRAPH = 'ms-dn_ds_graph'
ID_APPLY_BUTTOM = 'ms-apply-buttom'

def get_tab_mutation_stats(queries: Queries, data_source: DataSource):
    return dbc.CardBody(
            children=[
                get_samples_tab_left_bar(queries, data_source),
                html.Div(
                    className="one columns",
                    children=[html.Br()]),
                get_mutation_stats_tab_graphs()
        ])


def get_mutation_stats_tab_graphs():
    return html.Div(
        className="nine columns",
        children=[
            html.Br(),
            html.Div(children=[
                html.Div(id=ID_VARIANTS_PER_SAMPLE_GRAPH, className="six columns"),
                html.Div(id=ID_SUBSTITUTIONS_GRAPH, className="five columns"),
            ]),
            html.Div(children=[
                html.Div(id=ID_INDEL_LENGTH_GRAPH, className="six columns"),
                html.Div(id=ID_ANNOTATIONS_GRAPH, className="five columns"),
            ]),
        ])


def get_samples_tab_left_bar(queries: Queries, data_source: DataSource):
    return html.Div(
        className="two columns",
        children=[
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
            dcc.Markdown("""Select one or more genes"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_GENE,
                options=[{'label': g.name, 'value': g.name} for g in queries.get_genes()],
                value=None,
                multi=True
            ),
            html.Br(),
            dcc.Markdown("""Select a variant type for mutations per sample and top mutations"""),
            dcc.Dropdown(
                id=ID_DROPDOWN_VARIANT_TYPE,
                options=[
                    {'label': VariantType.SNV.name, 'value': VariantType.SNV.name},
                    {'label': VariantType.MNV.name, 'value': VariantType.MNV.name},
                    {'label': VariantType.INSERTION.name, 'value': VariantType.INSERTION.name},
                    {'label': VariantType.DELETION.name, 'value': VariantType.DELETION.name}
                ],
                value=None,
                multi=True
            ),
            html.Br(),
            html.Button('Apply', id=ID_APPLY_BUTTOM),
        ])


def set_callbacks_mutation_stats_tab(app, session: Session):

    queries = Queries(session=session)
    figures = MutationStatsFigures(queries)

    @app.callback(
        Output(ID_VARIANTS_PER_SAMPLE_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_GENE, 'value'),
            State(ID_DROPDOWN_VARIANT_TYPE, 'value')
        ],
        suppress_callback_exceptions=True
    )
    def update_variants_per_sample(_, data_source, genes, variant_types):
        return html.Div(children=figures.get_variants_per_sample_plot(
            data_source=data_source, genes=genes, variant_types=variant_types))

    @app.callback(
        Output(ID_SUBSTITUTIONS_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_GENE, 'value'),
            State(ID_DROPDOWN_VARIANT_TYPE, 'value')
        ],
        suppress_callback_exceptions=True
    )
    def update_substitutions(_, data_source, genes, variant_types):
        return html.Div(children=figures.get_substitutions_plot(
            data_source=data_source, genes=genes, variant_types=variant_types))

    @app.callback(
        Output(ID_INDEL_LENGTH_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_GENE, 'value')
        ],
        suppress_callback_exceptions=True
    )
    def update_indel_lengths(_, data_source, genes):
        return html.Div(children=figures.get_indels_lengths_plot(data_source=data_source, genes=genes))

    @app.callback(
        Output(ID_ANNOTATIONS_GRAPH, 'children'),
        inputs=[Input(ID_APPLY_BUTTOM, 'n_clicks')],
        state=[
            State(ID_DROPDOWN_DATA_SOURCE, 'value'),
            State(ID_DROPDOWN_GENE, 'value')
        ],
        suppress_callback_exceptions=True
    )
    def update_annotation(_, data_source, genes):
        return html.Div(children=figures.get_annotations_plot(data_source=data_source, genes=genes))