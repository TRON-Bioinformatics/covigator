import os
from datetime import datetime

import dash
import dash_core_components as dcc
import dash_html_components as html
from sqlalchemy import func, text
from sqlalchemy.orm import Session
from tenacity import wait_exponential, stop_after_attempt

from covigator import ENV_COVIGATOR_DASHBOARD_HOST, ENV_COVIGATOR_DASHBOARD_PORT, ENV_COVIGATOR_DASHBOARD_LOG_FILE
from covigator.dashboard.figures import get_accumulated_samples_by_country, get_variants_plot, get_circos_plot
from covigator.database.model import SampleEna, JobEna, JobStatus, Variant, VariantObservation, DataSource
from covigator.database.database import Database
import tenacity
import logzero
from logzero import logger

from covigator.database.queries import get_date_of_first_ena_sample, get_date_of_most_recent_ena_sample, \
    get_date_of_last_check, get_date_of_last_update

MISSING_VALUE = "-"


@tenacity.retry(wait=wait_exponential(multiplier=2, min=1, max=10), stop=stop_after_attempt(5))
def get_database():
    try:
        database = Database()
        session = database.get_database_session()
        stmt = text("SELECT 1")
        session.execute(stmt)
        logger.info("Database connected!")
    except Exception as e:
        logger.error("Connection to database failed, retrying...")
        raise e
    return database


def get_header():
    logo = "/assets/CoVigator_logo_txt_nobg.png"
    return html.Div(
            [
                html.Div(
                    [html.Img(src=logo, id="covigator-logo",
                              style={"height": "100px", "width": "auto", "margin-bottom": "25px",},)
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [html.Div([html.H1("Monitoring SARS-Cov-2 mutations")], style={"text-align": "center"})],
                    className="two column",
                    id="title",
                ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        )


def get_footer():
    return html.Footer(
        [
            html.Br(),
            html.Div(
                [
                    html.P(),
                    html.P("© 2021 TRON Mainz. All Rights Reserved"),
                    html.P()
                ],
                className="one-third column"
            ),
            html.Div(
                [
                    html.P(),
                    html.P([html.A("DATA PROTECTION", href="https://tron-mainz.de/data-protection/"), " | ",
                            html.A("IMPRINT", href="https://tron-mainz.de/imprint/")]),
                    html.P()
                ],
                className="one-third column"
            ),
            html.Br(),
        ],
        # this bit makes sure the footer sticks at the bottom
        style={"position": "relative", "bottom": "0", "width": "100%", "overflow": "hidden",
               "height": "100px"})


def get_tab_overview(session: Session):

    count_samples = session.query(JobEna).filter(JobEna.status == JobStatus.LOADED).count()
    count_countries = session.query(SampleEna).join(JobEna).filter(JobEna.status == JobStatus.LOADED).distinct(SampleEna.country)\
        .count()
    count_variants = session.query(Variant).count()
    count_variants_observed = session.query(VariantObservation).count()
    count_library_strategies = session.query(SampleEna.library_strategy, func.count(SampleEna.library_strategy)) \
        .join(JobEna).filter(JobEna.status == JobStatus.LOADED)\
        .group_by(SampleEna.library_strategy).all()
    count_instrument_model = session.query(SampleEna.instrument_model, func.count(SampleEna.instrument_model)) \
        .join(JobEna).filter(JobEna.status == JobStatus.LOADED) \
        .group_by(SampleEna.instrument_model).all()

    date_of_first_ena_sample = str(get_date_of_first_ena_sample(session))
    date_of_most_recent_ena_sample = str(get_date_of_most_recent_ena_sample(session))
    date_of_last_check_ena = str(get_date_of_last_check(session, data_source=DataSource.ENA))
    date_of_last_update_ena = str(get_date_of_last_update(session, data_source=DataSource.ENA))

    return dcc.Tab(label="About",
                        children=[html.Div(
                            [
                                html.H3("Covigator - Monitoring SARS-Cov-2 mutations"),
                                html.P("Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, "
                                       "necessitating preventive or therapeutic strategies and first steps towards an end to "
                                       "this pandemic were done with the approval of the first mRNA vaccines against SARS-CoV-2. "
                                       "We want to provide an interactive view on different types of variants that can be "
                                       "incorporated in global efforts to sustainably prevent or treat infections. Thus, we "
                                       "envision to help guiding global vaccine design efforts to overcome the threats of this "
                                       "pandemic."),
                                html.Br(),
                                html.P("If you want to cite us:"),
                                html.P([
                                    "Schrörs, B., Gudimella, R., Bukur, T., Rösler, T., Löwer, M., & Sahin, U. (2021). "
                                    "Large-scale analysis of SARS-CoV-2 spike-glycoprotein mutants demonstrates the "
                                    "need for continuous screening of virus isolates. BioRxiv, 2021.02.04.429765. ",
                                    html.A("https://doi.org/10.1101/2021.02.04.429765",
                                           href="https://doi.org/10.1101/2021.02.04.429765")],
                                    style={"font-style": "italic"}),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.H6("No. of Samples"), html.H6(count_samples)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("No. of countries"), html.H6(count_countries)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("No. of unique variants"), html.H6(count_variants)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("No. of variant observations"), html.H6(count_variants_observed)],
                                            className="mini_container",
                                        ),
                                    ],
                                    id="info-container-1",
                                    className="row container-display",
                                ),
                                html.Div(html.H4("European Nucleotide Archive (ENA)")),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.H6("First sample"),
                                             html.H6(_print_date(date_of_first_ena_sample))],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Most recent sample"),
                                             html.H6(_print_date(date_of_most_recent_ena_sample))],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Last checked"),
                                             html.H6(_print_date(date_of_last_check_ena))],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Last updated"),
                                             html.H6(_print_date(date_of_last_update_ena))],
                                            className="mini_container",
                                        ),
                                    ],
                                    id="info-container-2",
                                    className="row container-display",
                                ),
                                html.Div(html.H4("GISAID")),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.H6("First sample"), html.H6(MISSING_VALUE)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Most recent sample"), html.H6(MISSING_VALUE)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Last checked"), html.H6(MISSING_VALUE)],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Last updated"), html.H6(MISSING_VALUE)],
                                            className="mini_container",
                                        ),
                                    ],
                                    id="info-container-3",
                                    className="row container-display",
                                ),
                                html.Div(
                                    [
                                        html.Div(
                                            [html.H6("Sequencing technologies")] +
                                            [html.P("{} - {}".format(l, c)) for l, c in count_library_strategies],
                                            className="mini_container",
                                        ),
                                        html.Div(
                                            [html.H6("Sequencing instruments")] +
                                            [html.P("{} - {}".format(l, c)) for l, c in count_instrument_model],
                                            className="mini_container",
                                        ),
                                    ],
                                    id="info-container-4",
                                    className="row container-display",
                                ),
                                html.Br(),
                                html.Br(),
                                html.Br(),
                            ],
                            style={"text-align": "left"}
                        )])


def _print_date(date: datetime.date):
    return str(date) if date else MISSING_VALUE


def get_tab_samples(session: Session):
    figure = get_accumulated_samples_by_country(session)
    return dcc.Tab(label="Samples",
                        children=[
                            dcc.Graph(figure=figure)
                        ])


def get_tab_variants(session: Session):

    # TODO: parametrise the gene name
    gene_name = "S"
    needle_plot = get_variants_plot(session, gene_name=gene_name)
    circos_plot = get_circos_plot()

    return dcc.Tab(label="Variants",
                   children=[html.Div(
                       id='needleplot-body', className="row container-display",
                       children=[
                           html.Div(children=[
                               dcc.Markdown("""
                                 **Select a gene**
                                 """),
                               dcc.Dropdown(
                                   id='dropdown',
                                   options=[{'label': c, 'value': c} for c in ["S", "N", "E"]],
                                   value=gene_name,
                                   multi=False
                               ),
                           ], className="three columns"),
                           html.Div(children=needle_plot, className="nine columns")
                           ]
                   ),
                       circos_plot
                   ]
                   )


def serve_layout():

    session = database.get_database_session()
    header = get_header()
    footer = get_footer()
    tab_overview = get_tab_overview(session)
    tab_samples = get_tab_samples(session)
    tab_variants = get_tab_variants(session)

    # assemble tabs in dcc.Tabs object
    tabs = dcc.Tabs(children=[tab_overview, tab_samples, tab_variants])

    # create layout
    layout = html.Div(children=[header, tabs, footer])

    session.close()

    return layout


class CovigatorDashBoardInitialisationError(Exception):
    pass


def main(debug=False):
    log_file = os.getenv(ENV_COVIGATOR_DASHBOARD_LOG_FILE)
    if log_file is not None:
        logzero.logfile(log_file, maxBytes=1e6, backupCount=3)
    host = os.getenv(ENV_COVIGATOR_DASHBOARD_HOST, "0.0.0.0")
    try:
        port = int(os.getenv(ENV_COVIGATOR_DASHBOARD_PORT, "8050"))
    except ValueError as e:
        logger.exception(e)
        raise CovigatorDashBoardInitialisationError(e)
    app.run_server(debug=debug, host=host, port=port)


# creates the Dash application
app = dash.Dash(
    __name__,
    title="CoVigator",
    meta_tags=[
        # A description of the app, used by e.g.
        # search engines when displaying search results.
        {
            'name': 'description',
            'content': 'CoVigator - monitoring Sars-Cov-2 mutations'
        },
        # A tag that tells Internet Explorer (IE)
        # to use the latest renderer version available
        # to that browser (e.g. Edge)
        {
            'http-equiv': 'X-UA-Compatible',
            'content': 'IE=edge'
        },
        # A tag that tells the browser not to scale
        # desktop widths to fit mobile screens.
        # Sets the width of the viewport (browser)
        # to the width of the device, and the zoom level
        # (initial scale) to 1.
        #
        # Necessary for "true" mobile support.
        {
          'name': 'viewport',
          'content': 'width=device-width, initial-scale=1.0'
        }
    ]
)
database = get_database()
app.layout = serve_layout


if __name__ == '__main__':
    main(debug=True)
