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
from covigator.dashboard.tabs.dataset_gisaid import get_tab_dataset_gisaid
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

TAB_STYLE = {"color": "#003c78", 'margin-right': '15px'}

ID_TAB_CONTENT = "tab-content"
DOWNLOAD_TAB_ID = "download"
HELP_TAB_ID = "help"
INTRAHOST_MUTATIONS_TAB_ID = "subclonal-variants"
RECURRENT_MUTATIONS_TAB_ID = "variants"
SAMPLES_TAB_ID = "samples"
LINEAGES_TAB_ID = "lineages"
MUTATIONS_TAB_ID = "mutation-stats"
GISAID_DATASET_TAB_ID = "gisaid-dataset"
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
                                            children=html.A(html.Img(src="/assets/CoVigator_logo_txt_reg_no_bg.png",
                                                                     height="25px"), href="/"),
                                            className="ml-2",
                                            id="logo"
                                        ),
                                        #dbc.Col(html.Br(), className="ml-2"),
                                        #dbc.Col(html.Br(), className="ml-2"),
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
                                                    "Home", href="/", class_name="m-1",
                                                    style={'font-size' : '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "GISAID dataset", href="/gisaid",
                                                    style={'font-size' : '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "ENA dataset", href="/ena",
                                                    style={'font-size' : '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "Documentation", href="https://covigator.readthedocs.io/en/latest",
                                                    target="_blank",
                                                    style={'font-size': '150%', "color": "#003c78"}),
                                                dbc.DropdownMenuItem(
                                                    "Acknowledgements", href="/acknowledgements",
                                                    style={'font-size' : '150%', "color": "#003c78"}),
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


def set_callbacks(app, session: Session, content_folder):

    queries = Queries(session=session)

    MAIN_PAGE = "main"
    ENA_PAGE = DataSource.ENA
    GISAID_PAGE = DataSource.GISAID
    ACKNOWLEDGEMENTS_PAGE = "acknowledgements"

    def _get_page(url):
        if url in ["", "/"]:
            return MAIN_PAGE
        elif url == "/gisaid":
            return GISAID_PAGE
        elif url == "/ena":
            return ENA_PAGE
        elif url == "/acknowledgements":
            return ACKNOWLEDGEMENTS_PAGE
        else:
            raise ValueError("This URL does not exist")

    @app.callback(
        Output('tabs', "children"),
        Output('tabs', "active_tab"),
        [Input("url", "pathname")])
    def switch_page(url):
        page = _get_page(url)
        if page == MAIN_PAGE:
            # show overview with links
            return None, None
        elif page == GISAID_PAGE:
            # show gisaid tabs
            return [
                dbc.Tab(label="Overview", tab_id=GISAID_DATASET_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Lineages", tab_id=LINEAGES_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Mutation statistics", tab_id=MUTATIONS_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Recurrent mutations", tab_id=RECURRENT_MUTATIONS_TAB_ID, label_style=TAB_STYLE),
                ], GISAID_DATASET_TAB_ID
        elif page == ENA_PAGE:
            # show ena tabs
            return [
               dbc.Tab(label="Overview", tab_id=ENA_DATASET_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Lineages", tab_id=LINEAGES_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Mutation statistics", tab_id=MUTATIONS_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Recurrent mutations", tab_id=RECURRENT_MUTATIONS_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Intrahost mutations", tab_id=INTRAHOST_MUTATIONS_TAB_ID, label_style=TAB_STYLE),
                dbc.Tab(label="Download data", tab_id=DOWNLOAD_TAB_ID, label_style=TAB_STYLE)], ENA_DATASET_TAB_ID
        elif page == ACKNOWLEDGEMENTS_PAGE:
            # show ena tabs
            return [
                dbc.Tab(label="Acknowledgements", tab_id=HELP_TAB_ID, label_style={"color": "#003c78", 'display': 'none'})], HELP_TAB_ID

    @app.callback(
        Output('logo', "children"),
        [Input("url", "pathname")])
    def switch_logo(url):
        page = _get_page(url)
        if page == MAIN_PAGE or page == ACKNOWLEDGEMENTS_PAGE:
            return html.A(html.Img(src="/assets/CoVigator_logo_txt_reg_no_bg.png", height="80px"), href="/")
        elif page == GISAID_PAGE:
            return html.A(html.Img(src="/assets/CoVigator_logo_GISAID.png", height="80px"), href="/")
        elif page == ENA_PAGE:
            return html.A(html.Img(src="/assets/CoVigator_logo_ENA.png", height="80px"), href="/")

    @app.callback(
        Output('top-right-logo', "children"),
        [Input("url", "pathname")])
    def switch_lat_update(url):
        page = _get_page(url)
        if page == MAIN_PAGE:
            return None
        elif page == GISAID_PAGE:
            return dbc.Button(
                "last updated {date}".format(date=queries.get_last_update(DataSource.GISAID)),
                outline=True, color="dark", className="me-1",
                # 'background-color': '#b71300',
                style={"margin-right": "15px", 'font-size': '85%'})
        elif page == ENA_PAGE:
            return dbc.Button(
                "last updated {date}".format(date=queries.get_last_update(DataSource.ENA)),
                outline=True, color="dark", className="me-1",
                style={"margin-right": "15px", 'font-size': '85%'})

    @app.callback(
        Output(ID_TAB_CONTENT, "children"),
        [Input("tabs", "active_tab"), Input("url", "pathname")])
    def switch_tab(at, url):
        page = _get_page(url)
        try:
            if page == MAIN_PAGE:
                return get_tab_overview()
            elif at == ENA_DATASET_TAB_ID:
                return get_tab_dataset_ena(queries=queries)
            elif at == GISAID_DATASET_TAB_ID:
                return get_tab_dataset_gisaid(queries=queries)
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


def main(debug=False):
    Dashboard(config=Configuration()).start_dashboard(debug=debug)


if __name__ == '__main__':
    main(debug=True)
