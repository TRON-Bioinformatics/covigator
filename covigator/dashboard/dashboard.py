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
from covigator.dashboard.tabs.footer import get_footer
from covigator.dashboard.tabs.overview import get_tab_overview
from covigator.dashboard.tabs.samples import get_tab_samples, set_callbacks_samples_tab
from covigator.dashboard.tabs.variants import get_tab_variants, set_callbacks_variants_tab
from covigator.database.database import get_database
from logzero import logger
from covigator.database.queries import Queries

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
        layout = dbc.Card([
            dbc.CardHeader(
                dbc.Tabs([
                    dbc.Tab(label="Overview", tab_id=OVERVIEW_TAB_ID),
                    dbc.Tab(label="ENA dataset", tab_id=ENA_DATASET_TAB_ID),
                    dbc.Tab(label="GISAID dataset", tab_id=GISAID_DATASET_TAB_ID),
                    dbc.Tab(label="Samples", tab_id=SAMPLES_TAB_ID),
                    dbc.Tab(label="Recurrent variants", tab_id=VARIANTS_TAB_ID)],
                    id="tabs",
                    active_tab="overview",
                    card=True,
                )),
            dbc.CardBody(dcc.Loading(id="loading-1", children=[html.Div(id="content")], type="default")),
            dbc.CardFooter(footer)
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
        set_callbacks(app=app, session=session)
        set_callbacks_variants_tab(app=app, session=session)
        set_callbacks_samples_tab(app=app, session=session)
        return app


def set_callbacks(app, session: Session):

    queries = Queries(session=session)

    @app.callback(Output("content", "children"), [Input("tabs", "active_tab")])
    def switch_tab(at):
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
        return html.P("This shouldn't ever be displayed...")


def main(debug=False):
    Dashboard(config=Configuration()).start_dashboard(debug=debug)


if __name__ == '__main__':
    main(debug=True)
