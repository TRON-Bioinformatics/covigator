import dash_bootstrap_components as dbc
import pandas as pd
from dash import dcc
from dash import Input, Output, State, html
from covigator.database.queries import Queries
from covigator.dashboard.tabs import COLOR_STATUS, NEWS_PATTERN


def get_tab_overview(queries: Queries):
    present_icon = html.I(className="fas fa-check", style={"color": "green"})
    absent_icon = html.I(className="fas fa-times", style={"color": "red"})
    features = [
        (present_icon, present_icon),  # "Mutations per sample"
        (present_icon, present_icon),  # "Most frequent base substitutions"
        (present_icon, present_icon),  # "Indel length distribution"
        (present_icon, present_icon),  # "Most frequent mutation effects"
        (present_icon, present_icon),  # "Samples accumulation"
        (present_icon, present_icon),  # "Evolutionary pressure (dN/dS)"
        (present_icon, present_icon),  # "Instantaneous lineage prevalence"
        (present_icon, present_icon),  # "Top recurrent mutations"
        (present_icon, present_icon),  # "Genome/gene view"
        (present_icon, present_icon),  # "Co-occurrence clustering"
        (present_icon, absent_icon)  # "Top intrahost mutations"
    ]
    features = pd.DataFrame.from_records(features,
                                         columns=["ENA", "COVID-19"],
                                         index=["Mutations per sample", "Most frequent base substitutions",
                                                "Indel length distribution", "Most frequent mutation effects",
                                                "Samples accumulation", "Evolutionary pressure (dN/dS)",
                                                "Instantaneous lineage prevalence", "Top recurrent mutations",
                                                "Genome/gene view", "Co-occurrence clustering",
                                                "Top intrahost mutations"])
    features.index.set_names("Dashboard features", inplace=True)
    return dbc.CardBody([
            get_header(),
            dbc.Row([
                dbc.Col([None], width=2),
                dbc.Col([
                    html.Div(
                        children=[
                            dbc.Row([
                                dbc.Col([
                                    html.Br(),
                                    html.P(
                                    """
                                    CoVigator is a knowledge base and dashboard for SARS-CoV-2 mutations, built by integrating 
                                    a full variant calling pipeline. CoVigator dashboard provides readily available interpretation 
                                    of mutations along with respective lineages in an interactive visualization.
                                    """, style={"text-align": "justify"}),
                                    html.P(
                                        """
                                        The main goal for CoVigator is to provide a comprehensive resource with an up-to-date list 
                                        of mutations supporting the global efforts of finding emerging SARS-CoV-2 variants. CoVigator
                                        dashboard displays pre-calculated ranked recurrent and intrahost mutations through time, 
                                        which facilitates the monitoring of SARS-CoV-2 mutations.
                                        CoVigator process publicly available SARS-CoV-2 raw reads (i.e. FASTQs) from the 
                                        European Nucleotide Archive (ENA) and genome assemblies from the COVID-19 Data Portal. 
                                        CoVigator is open-data-friendly and allowing it to be adopted to other SARS-CoV-2 data sources. 
                                        """, style={"text-align": "justify"}),
                                    html.P([
                                        """
                                        CoVigator provides high-resolution SARS-CoV-2 mutations from genome assemblies and raw 
                                        reads allowing confirmation of evolutionary trends. CoVigator pipeline puts a special 
                                        emphasis on identification of SARS-CoV-2 intrahost mutations, thus reporting potential 
                                        SARS-CoV-2 variants of concern (VoC). CoVigator project is open sourced and made 
                                        available under the MIT license. The knowledge base and dashboard source is 
                                        hosted at """,
                                        html.A("github.",
                                               href="https://github.com/TRON-Bioinformatics/covigator",
                                               target="_blank"),
                                        """
                                        If you are interested in our work, please also have a read of our most recent 
                                        """,
                                        html.A("publication.",
                                               href="https://www.mdpi.com/1999-4915/15/6/1391",
                                               target="_blank")], style={"text-align": "justify"}),
                                    html.P(
                                        "Bukur, T., Riesgo-Ferreiro, P., Sorn, P., Gudimella, R., Hausmann, J., Rösler, T., "
                                        "Löwer, M., Schrörs, B., & Sahin, U. CoVigator — A Knowledge Base for Navigating SARS-CoV-2 Genomic Variants. "
                                        "Viruses. 2023; 15(6):1391.",
                                        style={"font-style": "italic", "margin-left": "50px"})
                                ], width=7),
                                dbc.Col(
                                    dbc.Card(
                                        [
                                            dbc.CardHeader([
                                                html.Div(
                                                    dbc.Row([
                                                        dbc.Col(html.H5("Covigator News"), width=6),
                                                        dbc.Col(
                                                            dbc.Button("see all covigator news",
                                                                       id="open",
                                                                       n_clicks=0,
                                                                       style={"justify": "end", "margin-right": "1%"}),
                                                            width=6),
                                                    ]),
                                                ),
                                                get_all_news(queries)
                                            ]),
                                            dbc.CardBody(
                                                html.Div(
                                                    children=get_news(queries, n=4)
                                                )
                                            )
                                        ], outline=True, style={"width": "40rem", "margin-left": "40px"},
                                    )
                                ),
                            ]),

                            dbc.CardBody(
                                dbc.Row([
                                    dbc.Col([html.Img(src="assets/wordcloud.png", style={"width": "100%"})]),
                                    dbc.Col([
                                        html.Br(),
                                        html.Div([
                                            dbc.Table.from_dataframe(features, striped=True, hover=True, index=True),
                                            #dbc.ListGroup(
                                            #    [
                                            #        dbc.ListGroupItem("Mutations per sample"),
                                            #        dbc.ListGroupItem("Most frequent base substitutions"),
                                            #        dbc.ListGroupItem("Indel length distribution"),
                                            #        dbc.ListGroupItem("Most frequent mutation effects"),
                                            #        dbc.ListGroupItem("Samples accumulation"),
                                            #        dbc.ListGroupItem("Evolutionary pressure (dN/dS)"),
                                            #        dbc.ListGroupItem("Instantaneous lineage prevalence"),
                                            #        dbc.ListGroupItem("Top recurrent mutations"),
                                            #        dbc.ListGroupItem("Genome/gene view"),
                                            #        dbc.ListGroupItem("Co-occurrence clustering"),
                                            #        dbc.ListGroupItem("Top intrahost mutations")
                                            #    ],
                                                #flush=True
                                            #)
                                        ]),
                                        ]),
                                    dbc.Col([
                                        html.Br(),
                                        dbc.Card(
                                            [
                                                dbc.CardImg(src="/assets/ENA_logo_2021.png", top=True,
                                                            style={"width": "23rem", "margin-left": "20px",
                                                                   "margin-top": "10px"}, ),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Button(
                                                            "Explore data derived from ENA", color="warning", href="/ena",
                                                            style={"margin-left": "20px", "margin-right": "20px",
                                                                   "font-size": 20}, ),
                                                    ]
                                                ),
                                            ],
                                            outline=False,
                                            style={"width": "40rem", "height": "15rem", "margin-left": "40px"},
                                        ),
                                        html.Br(),
                                        dbc.Card(
                                            [
                                                dbc.CardImg(src="/assets/CV19DP_logo_oneliner2.svg", top=True,
                                                            style={"width": "18rem", "margin-left": "20px",
                                                                   "margin-top": "10px"}, ),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Button(
                                                            "Explore data derived from the COVID-19 Data Portal sequences",
                                                            color="warning",
                                                            href="/covid19-portal",
                                                            style={"margin-left": "20px", "margin-right": "20px",
                                                                   "font-size": 20}, ),
                                                    ]
                                                ),
                                            ],
                                            outline=False,
                                            style={"width": "40rem", "height": "15rem", "margin-left": "40px"},
                                        )]
                                    ),
                                ])),
                            html.Br(),
                            html.Br(),
                            html.Br(),
                            html.Br(),
                        ],
                        style={"text-align": "left", "font-size": 16}),
                ], width=8),
                dbc.Col([None], width=2),
            ], align='center'),
        ])


def get_header():
    return dbc.Row([
                dbc.Col([None], width=2),
                dbc.Col([html.Div(
                    [html.Div([html.H1("Welcome to CoVigator Dashboard")], style={"text-align": "left"})],
                    id="title")], width=4),
                dbc.Col([None], width=2)
        ], id="header", className="row flex-display",)


def get_all_news(queries: Queries):
    news_items = get_news(queries, None)
    modal = dbc.Modal(
        [
            dbc.ModalHeader(dbc.ModalTitle("All Covigator News"), close_button=False),
            dbc.ModalBody(children=news_items),
            dbc.ModalFooter(
                dbc.Button(
                    "Close", id="close", className="ms-auto", n_clicks=0
                )
            ),
        ],
        id="modal",
        is_open=False,
        scrollable=True,
    )
    return modal


def get_news(queries: Queries, n=3):
    """
    Get news from database when website is (re)loaded and create children for news section
    cardboard.
    """
    newest_news = queries.get_top_news(n=n)
    children = []
    for this_news in newest_news.itertuples():
        children.append(
            dbc.Alert(children=[
                    dbc.Badge(children=[
                        this_news.published_date.strftime(NEWS_PATTERN),
                        " ",
                        this_news.message_type.name],
                        className="me-1"
                    ),
                    html.Br(),
                    dcc.Markdown(this_news.message_text)
                ],
                color=COLOR_STATUS[this_news.message_type.name]
            ))
    return children


def set_callbacks_news_section(app):
    @app.callback(
        Output("modal", "is_open"),
        [Input("open", "n_clicks"), Input("close", "n_clicks")],
        [State("modal", "is_open")],
    )
    def toggle_modal(n1, n2, is_open):
        if n1 or n2:
            return not is_open
        return is_open
