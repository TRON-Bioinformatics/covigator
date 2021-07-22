import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

from covigator.dashboard.tabs import get_mini_container, COLOR_OVERVIEW_MINI_CONTAINER, print_number
from covigator.database.queries import Queries


def get_tab_overview(queries: Queries):

    count_samples = queries.count_samples()
    count_countries = queries.count_countries()
    count_variants = queries.count_variants()

    return dbc.Card(
        dbc.CardBody([
            get_header(),
            html.Div(
               children=[
                   html.P(
                       "Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, "
                       "necessitating preventive or therapeutic strategies and first steps towards an end to "
                       "this pandemic were done with the approval of the first mRNA vaccines against SARS-CoV-2. "
                       "We want to provide an interactive view on different types of mutations that can be "
                       "incorporated in global efforts to sustainably prevent or treat infections. Thus, we "
                       "envision to help guiding global vaccine design efforts to overcome the threats of this "
                       "pandemic."),
                   dcc.Markdown("""
                       CoVigator loads publicly available SARS-CoV-2 DNA sequences from two systems: 

                       * European Nucleotide Archive (ENA) providing raw reads in FASTQ format
                       * Global Initiative on Sharing Avian Influenza Data (GISAID) providing assemblies in FASTA format
                       """),
                   html.Div(
                       className="row flex-display",
                       children=[
                       get_mini_container(
                           title="Samples",
                           value=print_number(count_samples),
                           color=COLOR_OVERVIEW_MINI_CONTAINER
                       ),
                       get_mini_container(
                           title="Countries",
                           value=print_number(count_countries),
                           color=COLOR_OVERVIEW_MINI_CONTAINER
                       ),
                       get_mini_container(
                           title="Mutations",
                           value=print_number(count_variants),
                           color=COLOR_OVERVIEW_MINI_CONTAINER
                       )
                   ]),
                   html.Br(),
                   html.P("If you want to cite us:"),
                   html.P([
                       "Schrörs, B., Gudimella, R., Bukur, T., Rösler, T., Löwer, M., & Sahin, U. (2021). "
                       "Large-scale analysis of SARS-CoV-2 spike-glycoprotein mutants demonstrates the "
                       "need for continuous screening of virus isolates. BioRxiv, 2021.02.04.429765. ",
                       ],
                       style={"font-style": "italic", "margin-left": "50px"}),
                   html.P(
                       html.A("https://doi.org/10.1101/2021.02.04.429765",
                              href="https://doi.org/10.1101/2021.02.04.429765"),
                       style={"text-indent": "50px"}
                   ),

                   html.Br(),
                   html.P(["The CoVigator analysis pipeline processes SARS-CoV-2 FASTQ or FASTA files into "
                           "annotated and normalized analysis ready VCF files. The pipeline is implemented "
                           "in the Nextflow framework (Di Tommaso, 2017). "
                           "It is open sourced under the MIT license, see the repository for "
                           "more details ",
                           html.A("https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline",
                                  href="https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline"),
                           "."]),
                   html.Br(),
                   html.P("If you would like to cite the pipeline:"),
                   html.P([
                       "Riesgo-Ferreiro, P., Sorn, P., & Bukur, T.. (2021, July 8). "
                       "TRON-Bioinformatics/covigator-ngs-pipeline: Release v0.5.0 (Version v0.5.0). "
                       "Zenodo. ",

                   ], style={"font-style": "italic", "margin-left": "50px"}),
                   html.P(
                       html.A("http://doi.org/10.5281/zenodo.5082444",
                              href="http://doi.org/10.5281/zenodo.5082444"),
                       style={"text-indent": "50px"}
                   ),
                   html.Br(),
               ],
               style={"text-align": "left", "font-size": 16}),
        ]))


def get_header():
    logo = "/assets/CoVigator_logo_txt_reg_no_bg.png"
    return html.Div(
            [
                html.Div(
                    [html.Img(src=logo, id="covigator-logo",
                              style={"height": "100px", "width": "auto", "margin-bottom": "10px"},)
                    ],
                ),
                html.Div(className="one column"),
                html.Div(
                    [html.Div([html.H1("Monitoring SARS-Cov-2 mutations")], style={"text-align": "left"})],
                    id="title",
                ),
            ],
            id="header",
            className="row flex-display",
        )
