# Run this app with `python dashboard.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
from sqlalchemy import func, text
from sqlalchemy.orm import Session
from tenacity import wait_exponential, stop_after_attempt

from covigator.dashboard.figures import get_accumulated_samples_by_country, get_variants_plot
from covigator.model import EnaRun, Job, JobStatus, Variant, VariantObservation
from covigator.database import Database
import tenacity
from logzero import logger


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

    count_samples = session.query(Job).filter(Job.status == JobStatus.LOADED).count()
    count_countries = session.query(EnaRun).join(Job).filter(Job.status == JobStatus.LOADED).distinct(EnaRun.country)\
        .count()
    count_variants = session.query(Variant).count()
    count_variants_observed = session.query(VariantObservation).count()
    count_library_strategies = session.query(EnaRun.library_strategy, func.count(EnaRun.library_strategy)) \
        .join(Job).filter(Job.status == JobStatus.LOADED)\
        .group_by(EnaRun.library_strategy).all()
    count_instrument_model = session.query(EnaRun.instrument_model, func.count(EnaRun.instrument_model)) \
        .join(Job).filter(Job.status == JobStatus.LOADED) \
        .group_by(EnaRun.instrument_model).all()

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
                                    id="info-container",
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
                                    id="info-container-2",
                                    className="row container-display",
                                ),
                                html.Br(),
                                html.Br(),
                                html.Br(),
                            ],
                            style={"text-align": "left"}
                        )])


def get_tab_samples(session: Session):
    figure = get_accumulated_samples_by_country(session)
    return dcc.Tab(label="Samples",
                        children=[
                            dcc.Graph(figure=figure)
                        ])


def get_tab_variants(session: Session):

    # TODO: parametrise the gene name
    gene_name = "S"
    figure = get_variants_plot(session, gene_name=gene_name)

    return dcc.Tab(label="Variants",
                   children=html.Div(
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
                           html.Div(children=figure, className="nine columns")

                           ]
                   )
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


def main(debug=False):
    app.run_server(debug=debug)


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
