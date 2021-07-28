import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from sqlalchemy.orm import Session
from dash.dependencies import Input, Output
import covigator
import covigator.configuration
from covigator.configuration import Configuration
from covigator.dashboard.tabs.dataset_ena import get_tab_dataset_ena
from covigator.dashboard.tabs.dataset_gisaid import get_tab_dataset_gisaid
from covigator.dashboard.tabs.download import set_callbacks_download_tab, get_tab_download
from covigator.dashboard.tabs.footer import get_footer
from covigator.dashboard.tabs.overview import get_tab_overview
from covigator.dashboard.tabs.samples import get_tab_samples, set_callbacks_samples_tab
from covigator.dashboard.tabs.subclonal_variants import get_tab_subclonal_variants, set_callbacks_subclonal_variants_tab
from covigator.dashboard.tabs.variants import get_tab_variants, set_callbacks_variants_tab
from covigator.database.database import get_database
from logzero import logger
from covigator.database.queries import Queries

DOWNLOAD_TAB_ID = "download"
SUBCLONAL_VARIANTS_TAB_ID = "subclonal-variants"
VARIANTS_TAB_ID = "variants"
SAMPLES_TAB_ID = "samples"
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
            dbc.Navbar([

                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.A(html.Img(
                            src="/assets/CoVigator_logo_txt_reg_no_bg.png", height="25px"),
                            href="https://covigator.tron-mainz.de/"), className="ml-2"
                        )
                    ],
                    align="center",
                    no_gutters=True,
                ),
                dbc.Row(
                    [
                        dbc.Col(html.A(html.Img(
                            src="/assets/tron_logo_without_text.png", height="22px"),
                            href="https://tron-mainz.de"), className="ml-2"),
                        dbc.Col(html.Br(), className="ml-2"),
                        dbc.Col(html.A(html.Img(
                            src="https://github.githubassets.com/images/modules/logos_page/Octocat.png", height="25px"),
                            href="https://github.com/TRON-bioinformatics/covigator"), className="ml-2"),
                        dbc.Col(html.Br(), className="ml-2"),

                    ],
                    align="center",
                    justify="end",
                    no_gutters=True,
                    style={'float': 'right', 'position': 'absolute', 'right': 0, 'text-align': 'right'}
                )
            ],
                color="white",
                dark=False,
            ),
            dbc.Card([
                dbc.CardHeader(
                    children=[
                        dbc.Tabs([
                            dbc.Tab(label="Overview", tab_id=OVERVIEW_TAB_ID),
                            dbc.Tab(label="ENA dataset", tab_id=ENA_DATASET_TAB_ID),
                            dbc.Tab(label="GISAID dataset", tab_id=GISAID_DATASET_TAB_ID),
                            dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID),
                            dbc.Tab(label="Recurrent clonal variants", tab_id=VARIANTS_TAB_ID),
                            dbc.Tab(label="Intrahost variants", tab_id=SUBCLONAL_VARIANTS_TAB_ID),
                            dbc.Tab(label="Download data", tab_id=DOWNLOAD_TAB_ID)],
                            id="tabs",
                            active_tab="overview",
                            card=True),

                    ]),
                dbc.CardBody(dcc.Loading(id="loading-1", children=[html.Div(id="content")], type="default")),
                dbc.CardFooter(footer)
            ])
        ])
        return layout

    def start_dashboard(self, debug=False):
        try:
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
        set_callbacks_subclonal_variants_tab(app=app, session=session)
        set_callbacks_download_tab(app=app, content_folder=self.config.content_folder)
        return app


def set_callbacks(app, session: Session, content_folder):

    queries = Queries(session=session)

    @app.callback(Output("content", "children"), [Input("tabs", "active_tab")])
    def switch_tab(at):
        logger.info("Changing tab...")
        try:
            if at == OVERVIEW_TAB_ID:
                return get_tab_overview(queries=queries)
            elif at == ENA_DATASET_TAB_ID:
                return get_tab_dataset_ena(queries=queries)
            elif at == GISAID_DATASET_TAB_ID:
                return get_tab_dataset_gisaid(queries=queries)
            elif at == SAMPLES_TAB_ID:
                return get_tab_samples(queries=queries)
            elif at == VARIANTS_TAB_ID:
                return get_tab_variants(queries=queries)
            elif at == SUBCLONAL_VARIANTS_TAB_ID:
                return get_tab_subclonal_variants(queries=queries)
            elif at == DOWNLOAD_TAB_ID:
                return get_tab_download(content_folder=content_folder)
            return html.P("This shouldn't ever be displayed...")
        except Exception as e:
            logger.exception(e)


def main(debug=False):
    Dashboard(config=Configuration()).start_dashboard(debug=debug)


if __name__ == '__main__':
    main(debug=True)
