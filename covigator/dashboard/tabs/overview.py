import dash_bootstrap_components as dbc
from dash import html
from sqlalchemy.orm import Session
from covigator.database.queries import Queries
from dash.dependencies import Output, Input, State

def get_tab_overview():

    return dbc.CardBody([
            get_header(),
            dbc.Row([
                dbc.Col([None], width=2),
                dbc.Col([
                    html.Div(
                        children=[
                            html.Br(),
                            html.P(
                                """
                                CoVigator is a knowledge base and dashboard for SARS-CoV-2 mutations, built by integrating 
                                a full variant calling pipeline. CoVigator dashboard provides readily available interpretation 
                                of mutations along with respective lineages in an interactive visualization.
                                """),
                            html.P(
                                """
                                The main goal for CoVigator is to provide a comprehensive resource with an up-to-date list 
                                of mutations supporting the global efforts of finding emerging SARS-CoV-2 variants. CoVigator
                                dashboard displays pre-calculated ranked recurrent and intrahost mutations through time, 
                                which facilitates the monitoring of SARS-CoV-2 mutations.
                                CoVigator process publicly available SARS-CoV-2 raw reads (i.e. FASTQs) from the 
                                European Nucleotide Archive (ENA) and genome assemblies from the COVID-19 Data Portal. 
                                CoVigator is open-data-friendly and allowing it to be adopted to other SARS-CoV-2 data sources. 
                                """),
                            html.Br(),
                            html.P("""
                               CoVigator provides high-resolution SARS-CoV-2 mutations from genome assemblies and raw 
                               reads allowing confirmation of evolutionary trends. CoVigator pipeline puts a special 
                               emphasis on identification of SARS-CoV-2 intrahost mutations, thus reporting potential 
                               SARS-CoV-2 variants of concern (VoC). CoVigator project is open sourced and made 
                               available under the MIT license. The knowledge base and dashboard source is 
                               hosted at https://github.com/TRON-Bioinformatics/covigator
                               """),
                            html.Br(),
                            html.P("""
                               If you are interested in our work, please also have a read of our most recent publication.
                               """),
                            html.P([
                                "Bukur, T., Riesgo-Ferreiro, P., Sorn, P., Gudimella, R., Hausmann, J., Rösler, T., "
                                "Löwer, M., Schrörs, B., & Sahin, U. CoVigator — A Knowledge Base for Navigating SARS-CoV-2 Genomic Variants. "
                                "Viruses. 2023; 15(6):1391.",
                            ], style={"font-style": "italic", "margin-left": "50px"}),


                            dbc.CardBody(
                                dbc.Row([
                                    dbc.Col([html.Img(src="assets/wordcloud.png", style={"width": "100%"})]),
                                    dbc.Col([
                                        html.Br(),
                                        html.Div([
                                            dbc.ListGroup(
                                                [
                                                    dbc.ListGroupItem("Mutations per sample"),
                                                    dbc.ListGroupItem("Most frequent base substitutions"),
                                                    dbc.ListGroupItem("Indel length distribution"),
                                                    dbc.ListGroupItem("Most frequent mutation effects"),
                                                    dbc.ListGroupItem("Samples accumulation"),
                                                    dbc.ListGroupItem("Evolutionary pressure (dN/dS)"),
                                                    dbc.ListGroupItem("Instantaneous lineage prevalence"),
                                                    dbc.ListGroupItem("Top recurrent mutations"),
                                                    dbc.ListGroupItem("Genome/gene view"),
                                                    dbc.ListGroupItem("Co-occurrence clustering"),
                                                    dbc.ListGroupItem("Top intrahost mutations")
                                                ],
                                                flush=True
                                            )
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
                                        ),
                                        html.Br(),
                                        dbc.Card(
                                            [
                                                dbc.CardBody(
                                                    dbc.Row([
                                                        dbc.Col([
                                                            html.Div(
                                                                children=[
                                                                    html.H2("Covigator News Section"),
                                                                    html.P(
                                                                        """
                                                                        Here you can find information about new data 
                                                                        releases, updates to Coivgator NGS pipeline or 
                                                                        Covigator python package
                                                                        """),
                                                                    html.Br(),
                                                                    html.Div(
                                                                        id="news_section",
                                                                        children=get_news(Session)
                                                                    )
                                                                ]
                                                            )
                                                        ])
                                                    ])
                                                )

                                            ],
                                            outline=False,
                                            style={"width": "40rem", "height": "15rem", "margin-left": "40px"},
                                        )
                                    ]),
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


def get_news(session: Session):
    queries = Queries(session=session)
    newest_news = queries.get_top_news(n=3)
    children = []
    #for this_news in newest_news.iterrows():
    #    news_item = []
    #    news_item.append(html.H3(this_news.message))
    #    news_item.append(html.P(this_news.publishing_date))
    #    children.extend(news_item)
    return children
