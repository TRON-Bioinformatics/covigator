import dash_bootstrap_components as dbc
from dash import html


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
                                Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, necessitating preventive or 
                                therapeutic strategies and first steps towards an end to this pandemic were done with the approval of the first mRNA 
                                vaccines against SARS-CoV-2. 
                                The accumulation of virus samples that have been sequenced in a short time frame is unprecedented.
                                This is the first pandemic recorded at a molecular level with such level of detail giving us the opportunity to develop
                                new tools for the monitoring of its evolution.
                                """),
                            html.P(
                                """
                                We want to provide an up-to-date interactive view on SARS-CoV-2 mutations to support global efforts in preventing or 
                                treating infections. 
                                Monitoring the appearance of relevant new mutations is key to enable a fast reaction to new strains and for that 
                                purpose we enable the exploration of these mutations and their annotations. 
                                Thus, we envision to help guiding global vaccine design efforts to overcome the threats of this pandemic.
                                """),
                            html.P(
                                """
                                CoVigator is a monitoring system for SARS-CoV-2 which integrates a full variant calling pipeline, 
                                a database that stores all relevant information about mutations in SARS-CoV-2 and finally a dashboard to enable 
                                visual analytics.
                                """),
                            html.Br(),
                            html.P("""
                       CoVigator loads publicly available SARS-CoV-2 DNA sequences from two systems: 
                       the European Nucleotide Archive (ENA) and the Global Initiative on Sharing Avian Influenza 
                       Data (GISAID). 
                       Some samples are present both in ENA and GISAID as some national initiatives are systematically
                       reporting to both databases. 
                       ENA provides the raw reads in FASTQ format and thus enables a higher resolution analysis into the 
                       SARS-CoV-2 mutations. Intrahost mutations are of particular interest. 
                       On the other hand, GISAID dataset is more extensive and represents best the geograhical and
                       temporal spread of SARS-CoV-2.
                       """),
                            html.Br(),
                            dbc.CardBody(
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card(
                                            [
                                                dbc.CardImg(src="/assets/gisaid.png", top=True,
                                                            style={"width": "18rem", "margin-left": "20px",
                                                                   "margin-top": "10px"}, ),
                                                dbc.CardBody(
                                                    [
                                                        dbc.Button(
                                                            "Explore data derived from GISAID", color="warning", href="/gisaid",
                                                            style={"margin-left": "20px", "margin-right": "20px",
                                                                   "font-size": 20}, ),
                                                    ]
                                                ),
                                            ],
                                            outline=False,
                                            style={"width": "40rem", "height": "15rem"},
                                        )
                                    ),
                                dbc.Col(
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
                                        style={"width": "40rem", "height": "15rem"},
                                    )
                                )])),
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
                    [html.Div([html.H1("Monitoring SARS-CoV-2 mutations")], style={"text-align": "left"})],
                    id="title")], width=5),
                dbc.Col([None], width=2)
        ], id="header", className="row flex-display",)