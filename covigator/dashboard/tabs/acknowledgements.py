import covigator
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html


def get_tab_acknowledgements():

    return dbc.CardBody([
        dbc.Row([
            dbc.Col([None], width=2),
            dbc.Col([
                html.H1("Acknowledgements"),
                dcc.Markdown("""
                **The CoVigator project** was developed at the Biomarker Development Center at TRON (Translational 
                Oncology at the University Medical Center of the Johannes Gutenberg University) gGmbH. 
                The project was kindly supported by Intel´s Pandemic Response Technology Initiative.
                """),
                html.Br(),
                html.P("If you want to cite us:"),
                html.P([
                    "Schrörs, B., Riesgo-Ferreiro, P., Sorn, P., Gudimella, R., Bukur, T., Rösler, T., "
                    "Löwer, M., & Sahin, U. (2021). Large-scale analysis of SARS-CoV-2 spike-glycoprotein mutants "
                    "demonstrates the need for continuous screening of virus isolates. PLOS ONE, 16(9), e0249254.",
                ],
                    style={"font-style": "italic", "margin-left": "50px"}),
                html.P(
                    html.A("https://doi.org/10.1101/2021.02.04.429765",
                           href="https://doi.org/10.1101/2021.02.04.429765", target="_blank"),
                    style={"text-indent": "50px"}
                ),
                html.Br(),
                html.P("If you would like to cite the pipeline:"),
                html.P([
                    "Riesgo-Ferreiro, P., Sorn, P., & Bukur, T.. (2021, July 8). "
                    "TRON-Bioinformatics/covigator-ngs-pipeline. Zenodo. ",

                ], style={"font-style": "italic", "margin-left": "50px"}),
                html.P(
                    html.A("https://doi.org/10.5281/zenodo.4906281",
                           href="https://doi.org/10.5281/zenodo.4906281", target="_blank"),
                    style={"text-indent": "50px"}
                ),
                html.Br(),

                html.P(html.B("GISAID Data")),
                html.P([
                    """We gratefully acknowledge all data contributors, i.e. the Authors and their Originating laboratories 
                    responsible for obtaining the specimens, and their Submitting laboratories for generating the 
                    genetic sequence and metadata and sharing via the GISAID Initiative""",
                    html.Sup("1"),
                    " on which this research is based."]),
                html.P([
                    html.Sup("1"),
                    """ Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: 
                    GISAID’s innovative contribution to global health. Global Challenges, 1:33-46. 
                    DOI: 10.1002/gch2.1018, PMCID: 31565258"""]),
                html.Br(),

                html.P(html.B("ENA Data")),
                html.P([
                    """We gratefully acknowledge the European Nucleotide Archive""",
                    html.Sup("2"),
                    " and all data contributors for sharing the raw reads on which this research is based."]),
                html.P([
                    html.Sup("2"),
                    """ Leinonen, R., Akhtar, R., Birney, E., Bower, L., Cerdeno-Tárraga, A., Cheng, Y., 
                    Cleland, I., Faruque, N., Goodgame, N., Gibson, R., Hoad, G., Jang, M., Pakseresht, N., 
                    Plaister, S., Radhakrishnan, R., Reddy, K., Sobhany, S., Hoopen, P. Ten, Vaughan, R., 
                    Zalunin V., Cochrane, G. (2011). The European nucleotide archive. Nucleic Acids Research, 
                    39(SUPPL. 1), D28. https://doi.org/10.1093/nar/gkq967"""]),
                html.Br(),

                html.P(html.B("About")),
                html.P("CoVigator dashboard version {}; analysis pipeline version {}".format(
                    covigator.VERSION, covigator.ANALYSIS_PIPELINE_VERSION
                )),
                html.P([
                    "Read the CoVigator documentation here ",
                    html.A("https://covigator.readthedocs.io",
                           href="https://covigator.readthedocs.io", target="_blank")
                ]),
                html.Br(),
                html.P("""
                       The CoVigator project is open sourced and made available under the MIT license.
                       We welcome any bug report, feature request or contribution through our GitHub repositories."""),
                html.P([
                    "The knowledge base and dashboard source is hosted at ",
                    html.A("https://github.com/TRON-Bioinformatics/covigator",
                           href="https://github.com/TRON-Bioinformatics/covigator", target="_blank")
                ]),
                html.P([
                    """
                    The analysis pipeline processes SARS-CoV-2 FASTQ or FASTA files into
                    annotated and normalized analysis ready VCF files. The pipeline is implemented
                    in the Nextflow framework (Di Tommaso, 2017), it is usable as an independent component
                    and we provide support to any user and to other viruses. The repository is hosted at
                    """,
                    html.A("https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline",
                           href="https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline", target="_blank")
                ]),

                html.Br(),
                html.P("We value your feedback!"),
                html.P("If you found an error or have a feature request, please, report it:"),
                html.Iframe(srcDoc="""
                        <a class="github-button" href="https://github.com/tron-bioinformatics/covigator/issues/new?assignees=&labels=bug&template=bug_report.md&title=" data-color-scheme="no-preference: light; light: light; dark: light;" data-size="large" data-show-count="false" aria-label="Issue tron-bioinformatics/covigator on GitHub">Issue</a>
                        <a class="github-button" href="https://github.com/tron-bioinformatics/covigator/issues/new?assignees=&labels=bug&template=feature_request.md&title=" data-color-scheme="no-preference: light; light: light; dark: light;" data-size="large" data-show-count="false" aria-label="Issue tron-bioinformatics/covigator on GitHub">Feature request</a>
                        <script async defer src="https://buttons.github.io/buttons.js"></script>
                       """, style={'border': 0, 'height': '50px'}),
                html.P("If you liked our work support us in GitHub:"),
                html.Iframe(srcDoc="""
                            <a class="github-button" href="https://github.com/tron-bioinformatics/covigator" data-color-scheme="no-preference: light; light: light; dark: light;" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star tron-bioinformatics/covigator on GitHub">Star</a>
                            <script async defer src="https://buttons.github.io/buttons.js"></script>
                           """, style={'border': 0, 'height': '50px'}),
                html.Br(),
                html.P("Follow #CoVigator in Twitter:"),
                html.A(html.Img(src="/assets/2021_twitter_logo_blue.png", height="25px"),
                       href="https://twitter.com/hashtag/CoVigator", target="_blank")

            ], width=8, style={"text-align": "left", "font-size": 16}),
            dbc.Col([None], width=2)]
            ),
    ])
