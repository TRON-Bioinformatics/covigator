import logging

import dash
from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
import logzero
from sqlalchemy.orm import Session
from dash.dependencies import Input, Output
import covigator
import covigator.configuration
from covigator.configuration import Configuration
from covigator.dashboard.tabs.acknowledgements import get_tab_acknowledgements
from covigator.dashboard.tabs.dataset_ena import get_tab_dataset_ena
from covigator.dashboard.tabs.dataset_covid19_portal import get_tab_dataset_covid19_portal
from covigator.dashboard.tabs.download import set_callbacks_download_tab, get_tab_download
from covigator.dashboard.tabs.footer import get_footer
from covigator.dashboard.tabs.lineages import set_callbacks_lineages_tab, get_tab_lineages
from covigator.dashboard.tabs.mutation_stats import get_tab_mutation_stats, set_callbacks_mutation_stats_tab
from covigator.dashboard.tabs.overview import get_tab_overview
from covigator.dashboard.tabs.samples import get_tab_samples, set_callbacks_samples_tab
from covigator.dashboard.tabs.intrahost_mutations import get_tab_subclonal_variants, set_callbacks_subclonal_variants_tab
from covigator.dashboard.tabs.recurrent_mutations import get_tab_variants, set_callbacks_variants_tab
from covigator.database.database import get_database
from logzero import logger

from covigator.database.model import DataSource
from covigator.database.queries import Queries

COVIGATOR_ENA_LOGO = "/assets/CoVigator_logo_ENA.png"
COVIGATOR_COVID19_LOGO = "/assets/CoVigator_logo_txt_reg_no_bg_covid19_portal.png"
COVIGATOR_LOGO = "/assets/CoVigator_logo_txt_reg_no_bg.png"
HOME_HREF = "/"
COVID_PORTAL_HREF = "/covid19-portal"
ENA_HREF = "/ena"
DOWNLOAD_HREF = "/download"
ACKNOWLEDGEMENTS_HREF = "/acknowledgements"
TAB_STYLE = {"color": "#003c78", 'margin-right': '15px'}

ID_TAB_CONTENT = "tab-content"
DOWNLOAD_TAB_ID = "download"
HELP_TAB_ID = "help"
INTRAHOST_MUTATIONS_TAB_ID = "subclonal-variants"
RECURRENT_MUTATIONS_TAB_ID = "variants"
SAMPLES_TAB_ID = "samples"
LINEAGES_TAB_ID = "lineages"
MUTATIONS_TAB_ID = "mutation-stats"
COVID19_PORTAL_DATASET_TAB_ID = "covid19-portal-dataset"
ENA_DATASET_TAB_ID = "ena-dataset"
OVERVIEW_TAB_ID = "overview"


class Dashboard:

    def __init__(self, config: Configuration, verbose=False):
        covigator.configuration.initialise_logs(config.logfile_dash)
        self.config = config
        # the connection to the database is created only once, but multiple sessions are used
        self.database = get_database(config=config, initialize=True, verbose=verbose)

    def serve_layout(self):
        logger.info("Serving layout")
        footer = get_footer()

        layout = html.Div(children=[
            dbc.Card([
                dbc.CardHeader(
                    children=[
                        dcc.Location(id='url', refresh=False),
                        dbc.Navbar(
                            [
                                dbc.Row([
                                        dbc.Col(
                                            children=html.A(html.Img(src=COVIGATOR_LOGO,
                                                                     height="80px"), href=HOME_HREF),
                                            className="ml-2",
                                            id="logo"
                                        ),
                                    ],
                                    align="left",
                                    className="g-0",
                                    style={'align': 'left'}
                                ),
                                dbc.Row(
                                    [
                                        dbc.Tabs(
                                            None,
                                            id="tabs",
                                            active_tab=SAMPLES_TAB_ID,
                                            style={'align': 'right', 'font-size': '140%'},
                                        )],
                                    align="right",
                                    style={'margin-left': '2%', }
                                ),
                                dbc.Row([
                                    dbc.Col(
                                        None,
                                        className="ml-2",
                                        id="top-right-logo",
                                        align="left",
                                        style={'margin-left': '2%', }
                                    ),
                                    dbc.Col(html.Br()),
                                    dbc.Col(
                                        dbc.DropdownMenu(
                                            label="Menu", children=[
                                                dbc.DropdownMenuItem(
                                                    "Home", href=HOME_HREF, class_name="m-1",
                                                    style={'font-size' : '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "ENA dashboard", href=ENA_HREF,
                                                    style={'font-size' : '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "COVID-19 Data Portal sequences dashboard", href=COVID_PORTAL_HREF,
                                                    style={'font-size': '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "Documentation", href="https://covigator.readthedocs.io/en/latest",
                                                    target="_blank",
                                                    style={'font-size': '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "Data download", href=DOWNLOAD_HREF,
                                                    style={'font-size': '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "Acknowledgements", href=ACKNOWLEDGEMENTS_HREF,
                                                    style={'font-size': '150%', "color": "#003c78"}),
                                            ],
                                            align_end=True,
                                            size="lg",
                                            toggle_style={
                                                "textTransform": "uppercase",
                                                "background": "#003c78",
                                                'font-size': '85%'
                                            },
                                            style={"margin-left": "15px"}
                                        ),
                                        style={}
                                    )],
                                    align="right",
                                    justify="end",
                                    style={'float': 'right', 'position': 'absolute', 'right': 0, 'text-align': 'right',
                                           },
                                    className="g-0 ms-auto flex-nowrap mt-3 mt-md-0",
                                )
                                ]
                        )],
                ),
                dbc.CardBody(dcc.Loading(id="loading-1", children=[html.Div(id=ID_TAB_CONTENT)], style={"height": "100%"})),
                dbc.CardFooter(footer)
            ])
        ])
        return layout

    def start_dashboard(self, debug=False):
        try:
            if debug:
                logzero.loglevel(level=logging.DEBUG)
            logger.info("Starting covigator dashboard")
            app = self.get_application()
            app.run_server(debug=debug, host=self.config.dash_host, port=self.config.dash_port)
        except Exception as e:
            logger.exception(e)
            raise e

    def get_application(self) -> dash.Dash:

        # creates the Dash application
        logger.info("Create the application")
        app = dash.Dash(
            name=__name__,
            title="CoVigator",
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
            meta_tags=[
                # A description of the app, used by e.g.
                # search engines when displaying search results.
                {
                    'name': 'description',
                    'content': 'CoVigator - monitoring SARS-CoV-2 mutations'
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
        # Warning pass the layout as a function, do not call it otherwise the application will serve a static dataset
        app.layout = self.serve_layout
        session = self.database.get_database_session()
        set_callbacks(app=app, session=session, content_folder=self.config.content_folder)
        set_callbacks_variants_tab(app=app, session=session)
        set_callbacks_samples_tab(app=app, session=session)
        set_callbacks_lineages_tab(app=app, session=session)
        set_callbacks_mutation_stats_tab(app=app, session=session)
        set_callbacks_subclonal_variants_tab(app=app, session=session)
        set_callbacks_download_tab(app=app, content_folder=self.config.content_folder)
        return app


MAIN_PAGE = "main"
ENA_PAGE = DataSource.ENA
COVID19_PORTAL_PAGE = DataSource.COVID19_PORTAL
ACKNOWLEDGEMENTS_PAGE = "acknowledgements"
DOWNLOAD_PAGE = "download"


def _get_page(url):
    if url in ["", HOME_HREF]:
        return MAIN_PAGE
    elif url == COVID_PORTAL_HREF:
        return COVID19_PORTAL_PAGE
    elif url == ENA_HREF:
        return ENA_PAGE
    elif url == ACKNOWLEDGEMENTS_HREF:
        return ACKNOWLEDGEMENTS_PAGE
    elif url == DOWNLOAD_HREF:
        return DOWNLOAD_PAGE
    else:
        raise ValueError("This URL does not exist")

def switch_page_callback(url):
    page = _get_page(url)
    if page == MAIN_PAGE:
        # show overview with links
        return None, None
    elif page == COVID19_PORTAL_PAGE:
        # show gisaid tabs
        return [
            dbc.Tab(label="Overview", tab_id=COVID19_PORTAL_DATASET_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Lineages", tab_id=LINEAGES_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Mutation statistics", tab_id=MUTATIONS_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Recurrent mutations", tab_id=RECURRENT_MUTATIONS_TAB_ID, label_style=TAB_STYLE),
        ], COVID19_PORTAL_DATASET_TAB_ID
    elif page == ENA_PAGE:
        # show ena tabs
        return [
            dbc.Tab(label="Overview", tab_id=ENA_DATASET_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Lineages", tab_id=LINEAGES_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Mutation statistics", tab_id=MUTATIONS_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Recurrent mutations", tab_id=RECURRENT_MUTATIONS_TAB_ID, label_style=TAB_STYLE),
            dbc.Tab(label="Intrahost mutations", tab_id=INTRAHOST_MUTATIONS_TAB_ID, label_style=TAB_STYLE)], \
            ENA_DATASET_TAB_ID
    elif page == ACKNOWLEDGEMENTS_PAGE:
        return [
            dbc.Tab(label="Acknowledgements", tab_id=HELP_TAB_ID,
                    label_style={"color": "#003c78", 'display': 'none'})], HELP_TAB_ID
    elif page == DOWNLOAD_PAGE:
        return [
            dbc.Tab(label="Download", tab_id=DOWNLOAD_TAB_ID,
                    label_style={"color": "#003c78", 'display': 'none'})], DOWNLOAD_TAB_ID


def switch_logo_callback(url):
    page = _get_page(url)
    if page in [MAIN_PAGE, ACKNOWLEDGEMENTS_PAGE, DOWNLOAD_PAGE]:
        return html.A(html.Img(src=COVIGATOR_LOGO, height="80px"), href="/")
    elif page == COVID19_PORTAL_PAGE:
        return html.A(html.Img(src=COVIGATOR_COVID19_LOGO, height="80px"), href="/")
    elif page == ENA_PAGE:
        return html.A(html.Img(src=COVIGATOR_ENA_LOGO, height="80px"), href="/")


def switch_lat_update_callback(url, queries):
    page = _get_page(url)
    if page == MAIN_PAGE:
        return None
    elif page == COVID19_PORTAL_PAGE:
        return dbc.Button(
            "last updated {date}".format(date=queries.get_last_update(DataSource.COVID19_PORTAL)),
            outline=True, color="dark", className="me-1",
            # 'background-color': '#b71300',
            style={"margin-right": "15px", 'font-size': '85%'})
    elif page == ENA_PAGE:
        return dbc.Button(
            "last updated {date}".format(date=queries.get_last_update(DataSource.ENA)),
            outline=True, color="dark", className="me-1",
            style={"margin-right": "15px", 'font-size': '85%'})


def switch_tab_callback(at, url, queries, content_folder):
    page = _get_page(url)
    try:
        if page == MAIN_PAGE:
            return get_tab_overview()
        elif at == ENA_DATASET_TAB_ID:
            return get_tab_dataset_ena(queries=queries)
        elif at == COVID19_PORTAL_DATASET_TAB_ID:
            return get_tab_dataset_covid19_portal(queries=queries)
        elif at == SAMPLES_TAB_ID:
            return get_tab_samples(queries=queries, data_source=page)
        elif at == MUTATIONS_TAB_ID:
            return get_tab_mutation_stats(queries=queries, data_source=page)
        elif at == RECURRENT_MUTATIONS_TAB_ID:
            return get_tab_variants(queries=queries, data_source=page)
        elif at == INTRAHOST_MUTATIONS_TAB_ID:
            return get_tab_subclonal_variants(queries=queries)
        elif at == DOWNLOAD_TAB_ID:
            return get_tab_download(content_folder=content_folder)
        elif at == HELP_TAB_ID:
            return get_tab_acknowledgements()
        elif at == LINEAGES_TAB_ID:
            return get_tab_lineages(queries=queries, data_source=page)
        return html.P("This shouldn't ever be displayed...")
    except Exception as e:
        logger.exception(e)


def set_callbacks(app, session: Session, content_folder):

    queries = Queries(session=session)

    @app.callback(
        Output('tabs', "children"),
        Output('tabs', "active_tab"),
        [Input("url", "pathname")])
    def switch_page(url):
        return switch_page_callback(url)

    @app.callback(
        Output('logo', "children"),
        [Input("url", "pathname")])
    def switch_logo(url):
        return switch_logo_callback(url)

    @app.callback(
        Output('top-right-logo', "children"),
        [Input("url", "pathname")])
    def switch_lat_update(url):
        return switch_lat_update_callback(url=url, queries=queries)

    @app.callback(
        Output(ID_TAB_CONTENT, "children"),
        [Input("tabs", "active_tab"), Input("url", "pathname")])
    def switch_tab(at, url):
        return switch_tab_callback(at=at, url=url, queries=queries, content_folder=content_folder)


def main(debug=False):
    Dashboard(config=Configuration()).start_dashboard(debug=debug)


if __name__ == '__main__':
    main(debug=True)
