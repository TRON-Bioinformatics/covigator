from datetime import datetime

import dash_core_components as dcc
import dash_html_components as html

from covigator.dashboard.tabs import TAB_STYLE, TAB_SELECTED_STYLE, MISSING_VALUE, get_mini_container
from covigator.database.model import DataSource
from covigator.database.queries import Queries


def get_tab_overview(queries: Queries):

    count_samples_ena = queries.count_samples(source=DataSource.ENA.name)
    count_samples_gisaid = queries.count_samples(source=DataSource.GISAID.name)
    count_countries = queries.count_countries()
    count_variants = queries.count_variants()
    count_insertions = queries.count_insertions()
    count_deletions = queries.count_deletions()
    count_variants_observed = queries.count_variant_observations()
    count_subclonal_variants_observed = queries.count_subclonal_variant_observations()
    # count_library_strategies = session.query(SampleEna.library_strategy, func.count(SampleEna.library_strategy)) \
    #    .join(JobEna).filter(JobEna.status == JobStatus.LOADED)\
    #    .group_by(SampleEna.library_strategy).all()
    # count_instrument_model = session.query(SampleEna.instrument_model, func.count(SampleEna.instrument_model)) \
    #    .join(JobEna).filter(JobEna.status == JobStatus.LOADED) \
    #    .group_by(SampleEna.instrument_model).all()

    date_of_first_ena_sample = queries.get_date_of_first_sample(source=DataSource.ENA)
    date_of_first_gisaid_sample = queries.get_date_of_first_sample(source=DataSource.GISAID)
    date_of_most_recent_ena_sample = queries.get_date_of_most_recent_sample(source=DataSource.ENA)
    date_of_most_recent_gisaid_sample = queries.get_date_of_most_recent_sample(source=DataSource.GISAID)
    date_of_last_update_ena = queries.get_date_of_last_update(data_source=DataSource.ENA)
    date_of_last_update_gisaid = queries.get_date_of_last_update(data_source=DataSource.GISAID)

    return dcc.Tab(label="About",
                   style=TAB_STYLE,
                   selected_style=TAB_SELECTED_STYLE,
                   children=[
                       html.Br(),
                       get_header(),
                       html.Div(
                           children=[
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
                                   children=[
                                       get_mini_container(
                                           title="No. of Samples",
                                           value=print_number(count_samples_ena + count_samples_gisaid)),
                                       get_mini_container(title="No. of countries", value=print_number(count_countries)),
                                       html.Div(
                                           children=[
                                               html.H6("No. of unique mutations"),
                                               html.H6(print_number(count_variants)),
                                               html.P("{} SNVs".format(
                                                   print_number(count_variants - count_insertions - count_deletions))),
                                               html.P("{} insertions".format(print_number(count_insertions))),
                                               html.P("{} deletions".format(print_number(count_deletions)))
                                           ],
                                           className="mini_container",
                                       ),
                                       get_mini_container(
                                           title="No. of mutation observations",
                                           value=print_number(count_variants_observed)),
                                       get_mini_container(
                                           title="No. of intrahost mutation observations",
                                           value=print_number(count_subclonal_variants_observed)),
                                   ],
                                   id="info-container-1",
                                   className="row container-display",
                               ),
                               html.Div(html.H4("European Nucleotide Archive (ENA)")),
                               html.Div(
                                   children=[
                                       get_mini_container(
                                           title="No. of Samples", value=print_number(count_samples_ena)),
                                       get_mini_container(
                                           title="First sample", value=print_date(date_of_first_ena_sample)),
                                       get_mini_container(
                                           title="Most recent sample", value=print_date(date_of_most_recent_ena_sample)),
                                       get_mini_container(
                                           title="Last updated", value=print_date(date_of_last_update_ena))
                                   ],
                                   id="info-container-2",
                                   className="row container-display",
                               ),
                               html.Div(html.H4("GISAID")),
                               html.Div(
                                   [
                                       get_mini_container(
                                           title="No. of Samples", value=print_number(count_samples_gisaid)),
                                       get_mini_container(
                                           title="First sample", value=print_date(date_of_first_gisaid_sample)),
                                       get_mini_container(
                                           title="Most recent sample",
                                           value=print_date(date_of_most_recent_gisaid_sample)),
                                       get_mini_container(
                                           title="Last updated", value=print_date(date_of_last_update_gisaid))
                                   ],
                                   id="info-container-3",
                                   className="row container-display",
                               ),
                               html.Div(
                                   [
                                       # html.Div(
                                       #    [html.H6("Sequencing technologies")] +
                                       #    [html.P("{} - {}".format(l, c)) for l, c in count_library_strategies],
                                       #    className="mini_container",
                                       # ),
                                       # html.Div(
                                       #    [html.H6("Sequencing instruments")] +
                                       #    [html.P("{} - {}".format(l, c)) for l, c in count_instrument_model],
                                       #    className="mini_container",
                                       # ),
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


def print_date(date: datetime.date):
    return str(date) if date is not None else MISSING_VALUE


def print_number(value):
    return '{:,}'.format(value)
