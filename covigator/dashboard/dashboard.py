from datetime import datetime, timedelta
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input
import dash_table
import covigator
import covigator.configuration
from covigator.configuration import Configuration
from covigator.database.model import DataSource
from covigator.database.database import session_scope, get_database
from logzero import logger
from covigator.database.queries import Queries
from covigator.dashboard.figures import Figures

MONTH_PATTERN = "%Y-%m"
MISSING_VALUE = "-"


class Dashboard:

    tab_style = {
        'borderBottom': '1px solid #d6d6d6',
        'padding': '6px',
        'fontWeight': 'bold'
    }

    tab_selected_style = {
        'borderTop': '1px solid #d6d6d6',
        'borderBottom': '1px solid #d6d6d6',
        'backgroundColor': '#119DFF',
        'color': 'white',
        'padding': '6px'
    }

    def __init__(self, config: Configuration, verbose=False):
        covigator.configuration.initialise_logs(config.logfile_dash)
        self.config = config
        # the connection to the database is created only once, but multiple sessions are used
        self.database = get_database(config=config, initialize=True, verbose=verbose)
        self.queries = None
        self.figures = None

    def get_header(self):
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

    def get_footer(self):
        return html.Footer(
            [
                html.Br(),
                html.Div(
                    [
                        html.P(),
                        html.P("Covigator {} © 2021 TRON Mainz. All Rights Reserved".format(covigator.VERSION)),
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

    def get_tab_overview(self):

        count_samples = self.queries.count_ena_samples()
        count_countries = self.queries.count_countries()
        count_variants = self.queries.count_variants()
        count_variants_observed = self.queries.count_variant_observations()
        #count_library_strategies = session.query(SampleEna.library_strategy, func.count(SampleEna.library_strategy)) \
        #    .join(JobEna).filter(JobEna.status == JobStatus.LOADED)\
        #    .group_by(SampleEna.library_strategy).all()
        #count_instrument_model = session.query(SampleEna.instrument_model, func.count(SampleEna.instrument_model)) \
        #    .join(JobEna).filter(JobEna.status == JobStatus.LOADED) \
        #    .group_by(SampleEna.instrument_model).all()

        date_of_first_ena_sample = self.queries.get_date_of_first_ena_sample()
        date_of_most_recent_ena_sample = self.queries.get_date_of_most_recent_ena_sample()
        date_of_last_check_ena = self.queries.get_date_of_last_check(data_source=DataSource.ENA)
        date_of_last_update_ena = self.queries.get_date_of_last_update(data_source=DataSource.ENA)

        return dcc.Tab(label="About", style=self.tab_style, selected_style=self.tab_selected_style,
                            children=[
                                html.Br(),
                                self.get_header(),
                                html.Div(
                                [
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
                                                 html.H6(self._print_date(date_of_first_ena_sample))],
                                                className="mini_container",
                                            ),
                                            html.Div(
                                                [html.H6("Most recent sample"),
                                                 html.H6(self._print_date(date_of_most_recent_ena_sample))],
                                                className="mini_container",
                                            ),
                                            html.Div(
                                                [html.H6("Last checked"),
                                                 html.H6(self._print_date(date_of_last_check_ena))],
                                                className="mini_container",
                                            ),
                                            html.Div(
                                                [html.H6("Last updated"),
                                                 html.H6(self._print_date(date_of_last_update_ena))],
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
                                            #html.Div(
                                            #    [html.H6("Sequencing technologies")] +
                                            #    [html.P("{} - {}".format(l, c)) for l, c in count_library_strategies],
                                            #    className="mini_container",
                                            #),
                                            #html.Div(
                                            #    [html.H6("Sequencing instruments")] +
                                            #    [html.P("{} - {}".format(l, c)) for l, c in count_instrument_model],
                                            #    className="mini_container",
                                            #),
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

    def _print_date(self, date: datetime.date):
        return str(date) if date is not None else MISSING_VALUE

    def get_tab_samples(self):
        figure = self.figures.get_accumulated_samples_by_country_plot()
        return dcc.Tab(label="Samples", style=self.tab_style, selected_style=self.tab_selected_style,
                            children=[
                                dcc.Graph(figure=figure)
                            ])

    def get_tab_variants(self):

        genes = self.queries.get_genes()
        months = self.queries.get_sample_months(MONTH_PATTERN)
        today = datetime.now()
        today_formatted = today.strftime(MONTH_PATTERN)
        oneyearago = today - timedelta(days=356)
        oneyearago_formatted = oneyearago.strftime(MONTH_PATTERN)

        return dcc.Tab(label="Variants", style=self.tab_style, selected_style=self.tab_selected_style,
                       children=[html.Div(
                           id='needleplot-body', className="row container-display",
                           style={'overflow': 'scroll'}, # 'top': 0, 'bottom': 0, position: fixed
                           children=[
                               html.Div(children=[
                                   html.Br(),
                                   dcc.Markdown("""
                                     Select a gene
                                     """),
                                   dcc.Dropdown(
                                       id='dropdown-gene',
                                       options=[{'label': c, 'value': c} for c in genes],
                                       value=None,
                                       multi=False
                                   ),
                                   html.Br(),
                                   dcc.Markdown("""
                                            **Top occurring variants**
                                             
                                             Number of top occurring variants
                                             """),
                                   dcc.Slider(
                                       id='slider-top-variants',
                                       min=10,
                                       max=100,
                                       step=10,
                                       value=10,
                                       dots=True,
                                       marks={i: '{}'.format(i) for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]},
                                       tooltip=dict(always_visible=False, placement="right")
                                   ),
                                   html.Br(),
                                   dcc.Markdown("""
                                                     Metric to measure abundance of variants per month
                                                     """),
                                   dcc.Dropdown(
                                       id='dropdown-top-variants-metric',
                                       options=[
                                           {'label': 'Count', 'value': 'count'},
                                           {'label': 'Frequency', 'value': 'frequency_by_month'}],
                                       value='count',
                                       clearable=False,
                                       multi=False
                                   ),
                                   html.Br(),
                                   dcc.Markdown("""
                                             Select a start and end date
                                             """),
                                   html.Div(children=[
                                       html.Div(
                                           dcc.Dropdown(
                                               id='dropdown-date-range-start',
                                               options=[{'label': c, 'value': c} for c in months],
                                               value=oneyearago_formatted,
                                               multi=False,
                                               clearable=False
                                           ), className="six column"),
                                       html.Div(id="dropdown-date-range-end-div",
                                                children=dcc.Dropdown(
                                                    id='dropdown-date-range-end',
                                                    options=[{'label': c, 'value': c} for c in months],
                                                    value=today_formatted,
                                                    multi=False,
                                                    clearable=False
                                                ),
                                                className="six column"),
                                   ], className="row container-display"),
                                   html.Br(),
                                   dcc.Markdown("""
                                                    **Genome view**
                                                     
                                                     Bin size
                                                     """),
                                   dcc.Slider(
                                       id='slider-bin-size',
                                       min=5,
                                       max=400,
                                       step=5,
                                       value=50,
                                       dots=False,
                                       marks={i: '{}'.format(i) for i in [10, 50, 100, 200, 300, 400]},
                                       tooltip=dict(always_visible=False, placement="right")
                                   ),
                                   html.Br(),
                                   dcc.Markdown("""
                                   **Co-occurrence matrix**
                                   
                                   Metric to assess paiwise co-occurrence
                                   """),
                                   dcc.Dropdown(
                                       id='dropdown-heatmap-metric',
                                       options=[{'label': "Count", 'value': "count"},
                                                {'label': "Frequency", 'value': "frequency"},
                                                {'label': "Jaccard index", 'value': "jaccard"},
                                                ],
                                       value="jaccard",
                                       clearable=False,
                                       multi=False
                                   ),
                                   dcc.Markdown("""
                                   Minimum number of pairwise co-occurrences
                                   """),
                                   dcc.Slider(
                                       id='slider-min-cooccurrences',
                                       min=1,
                                       max=100,
                                       step=5,
                                       value=20,
                                       dots=True,
                                       marks={i: '{}'.format(i) for i in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]},
                                       tooltip=dict(always_visible=False, placement="right")
                                   ),
                               ], className="two columns"),
                               html.Div(children=[
                                   html.Br(),
                                   html.Div(id='top-occurring-variants', children=dash_table.DataTable(id="top-occurring-variants-table")),
                                   html.Br(),
                                   html.Div(id='needle-plot'),
                                   html.Br(),
                                   html.Div(id='cooccurrence-heatmap'),
                                   html.Br(),
                                   html.Div(id='variants-mds'),
                                   html.Br(),
                               ], className="ten columns", style={'overflow': 'scroll', "height": "900px"},)
                           ]),
                       ]
                       )

    def serve_layout(self):
        logger.info("Serving layout")
        with session_scope(database=self.database) as session:
            self.queries = Queries(session=session)
            self.figures = Figures(self.queries)
            footer = self.get_footer()
            tab_overview = self.get_tab_overview()
            tab_samples = self.get_tab_samples()
            tab_variants = self.get_tab_variants()

            # assemble tabs in dcc.Tabs object
            tabs = dcc.Tabs(children=[tab_overview, tab_samples, tab_variants],
                            style={'height': '44px'})
                            #style={'margin': '0 0 0 0', 'padding-top': '2px'})

            # create layout
            layout = html.Div(children=[tabs, footer]) #, style={'margin-top': '0px', 'padding-top': '0px'})

        return layout

    def start_dashboard(self, debug=False):
        try:
            logger.info("Starting covigator dashboard")
            app = self.get_application()
            self.set_callbacks(app)
            app.run_server(debug=debug, host=self.config.dash_host, port=self.config.dash_port)
        except Exception as e:
            logger.exception(e)
            raise e

    def set_callbacks(self, app):
        @app.callback(
            Output('top-occurring-variants', 'children'),
            Input('slider-top-variants', 'value'),
            Input('dropdown-gene', 'value'),
            Input('dropdown-date-range-start', 'value'),
            Input('dropdown-date-range-end', 'value'),
            Input('dropdown-top-variants-metric', 'value')
        )
        def update_top_occurring_variants(top_variants, gene_name, date_range_start, date_range_end, metric):
            return html.Div(children=self.figures.get_top_occurring_variants_plot(
                top=top_variants, gene_name=gene_name, date_range_start=date_range_start,
                date_range_end=date_range_end, metric=metric))

        @app.callback(
            Output('needle-plot', 'children'),
            Input('dropdown-gene', 'value'),
            Input('top-occurring-variants-table', "derived_virtual_data"),
            Input('top-occurring-variants-table', "derived_virtual_selected_rows"),
            Input('slider-bin-size', 'value'),
        )
        def update_needle_plot(gene_name, rows, selected_rows_indices, bin_size):
            if gene_name is not None:
                selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
                plot = html.Div(
                    children=self.figures.get_variants_plot(
                        gene_name=gene_name, selected_variants=selected_rows, bin_size=bin_size))
            else:
                plot = html.Div(
                    children=self.figures.get_variants_abundance_plot(bin_size=bin_size))
            return plot

        @app.callback(
            Output('dropdown-date-range-end-div', 'children'),
            Input('dropdown-date-range-start', 'value'))
        def update_dropdown_end_date(start_date):
            today = datetime.now()
            today_formatted = today.strftime(MONTH_PATTERN)
            months = [m for m in self.queries.get_sample_months(MONTH_PATTERN) if m >= start_date]
            return dcc.Dropdown(
                id='dropdown-date-range-end',
                options=[{'label': c, 'value': c} for c in months],
                value=today_formatted,
                multi=False,
                clearable=False
            )

        @app.callback(
            Output('cooccurrence-heatmap', 'children'),
            Input('dropdown-gene', 'value'),
            Input('top-occurring-variants-table', "derived_virtual_data"),
            Input('top-occurring-variants-table', "derived_virtual_selected_rows"),
            Input('dropdown-heatmap-metric', 'value'),
            Input('slider-min-cooccurrences', 'value'),

        )
        def update_cooccurrence_heatmap(gene_name, rows, selected_rows_indices, metric, min_occurrences):

            selected_rows = [rows[s] for s in selected_rows_indices] if selected_rows_indices else None
            plot = html.Div(children=self.figures.get_cooccurrence_heatmap(
                gene_name=gene_name, selected_variants=selected_rows, metric=metric, min_occurrences=min_occurrences))
            return plot

        @app.callback(
            Output('variants-mds', 'children'),
            Input('dropdown-gene', 'value'),
        )
        def update_variants_mds(gene_name):
            # TODO pass selected variants
            plot = html.Div(children=self.figures.get_variants_clustering(gene_name=gene_name))
            return plot

    def get_application(self) -> dash.Dash:

        # creates the Dash application
        logger.info("Create the application")
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
        # Warning pass the layout as a function, do not call it otherwise the application will serve a static dataset
        app.layout = self.serve_layout
        return app


def main(debug=False):
    Dashboard(config=Configuration()).start_dashboard(debug=debug)


if __name__ == '__main__':
    main(debug=True)
