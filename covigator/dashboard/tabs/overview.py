import functools

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

from covigator.dashboard.tabs import get_mini_container, print_number
from covigator.database.queries import Queries


@functools.lru_cache()
def get_tab_overview(queries: Queries):

    count_samples = queries.count_samples()
    count_countries = queries.count_countries()
    count_variants = queries.count_variants()

    return dbc.CardBody([
            get_header(),
            html.Div(
               children=[
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
                   html.Div(
                       html.Span(
                           children=[
                               get_mini_container(
                                   title="Samples",
                                   value=print_number(count_samples)
                               ),
                               get_mini_container(
                                   title="Countries",
                                   value=print_number(count_countries)
                               ),
                               get_mini_container(
                                   title="Mutations",
                                   value=print_number(count_variants)
                               )
                           ])
                   ),
                   html.Br(),
                   dcc.Markdown("""
                       CoVigator loads publicly available SARS-CoV-2 DNA sequences from two systems: 

                       * [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena) providing raw reads in FASTQ format
                       * [Global Initiative on Sharing Avian Influenza Data (GISAID)](https://www.gisaid.org/) providing assemblies in FASTA format
                       """),
                   html.Br(),
                   html.P("""
                    There is certain overlap in the samples present in ENA and GISAID as some national initiatives are systematically 
                    reporting to both databases. ENA enables a higher resolution into the SARS-CoV-2 mutation details through the individual 
                    reads. This allows us to annotate mutations with a Variant Allele Frequency (VAF) and explore intrahost 
                    mutations. On the other hand, while we load all of the GISAID database in CoVigator, we only process the Illumina 
                    samples from ENA. This means excluding all of the Oxford Nanopore samples and hence having a partial view of all the 
                    available data.
                   """),
                   html.Br(),
                   html.P("If you want to cite us:"),
                   html.P([
                       "Schrörs, B., Riesgo-Ferreiro, P., Sorn, P., Gudimella, R., Bukur, T., Rösler, T., "
                       "Löwer, M., & Sahin, U. (2021). Large-scale analysis of SARS-CoV-2 spike-glycoprotein mutants "
                       "demonstrates the need for continuous screening of virus isolates. PLOS ONE, 16(9), e0249254. "
                       "https://doi.org/10.1371/journal.pone.0249254",
                       ],
                       style={"font-style": "italic", "margin-left": "50px"}),
                   html.P(
                       html.A("https://doi.org/10.1101/2021.02.04.429765",
                              href="https://doi.org/10.1101/2021.02.04.429765",  target="_blank"),
                       style={"text-indent": "50px"}
                   ),
                   html.Br(),
                   html.P("If you would like to cite the pipeline:"),
                   html.P([
                       "Riesgo-Ferreiro, P., Sorn, P., & Bukur, T.. (2021, July 8). "
                       "TRON-Bioinformatics/covigator-ngs-pipeline: Release v0.5.0 (Version v0.5.0). "
                       "Zenodo. ",

                   ], style={"font-style": "italic", "margin-left": "50px"}),
                   html.P(
                       html.A("http://doi.org/10.5281/zenodo.5082444",
                              href="http://doi.org/10.5281/zenodo.5082444",  target="_blank"),
                       style={"text-indent": "50px"}
                   ),
                   html.Br(),
               ],
               style={"text-align": "left", "font-size": 16}),
        ])


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
                    [html.Div([html.H1("Monitoring SARS-CoV-2 mutations")], style={"text-align": "left"})],
                    id="title",
                ),
            ],
            id="header",
            className="row flex-display",
        )
