import functools

import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html

import covigator
from covigator.dashboard.tabs import get_mini_container, print_number
from covigator.database.queries import Queries


def get_tab_help():

    return dbc.CardBody([
            html.Div(
               children=[
                   html.P(html.B("CoVigator version {}".format(covigator.VERSION))),
                   html.Br(),
                   dcc.Markdown("Read the CoVigator documentation here https://covigator.readthedocs.io"),
                   html.Br(),
                   dcc.Markdown("""
                       The CoVigator project is open sourced and made available under the MIT license.
                       We welcome any bug report, feature request or contribution through our GitHub repositories.

                       * The CoVigator knowledge base and dashboard are available at 
                       https://github.com/TRON-Bioinformatics/covigator.
                       * The CoVigator analysis pipeline processes SARS-CoV-2 FASTQ or FASTA files into
                       annotated and normalized analysis ready VCF files. The pipeline is implemented
                       in the Nextflow framework (Di Tommaso, 2017), it is usable as an independent component
                       and we provide support to any user. The repository is available at 
                       https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline
                       
                       The CoVigator project was developed at the Biomarker Development Center at TRON (Translational 
                       Oncology at the University Medical Center of the Johannes Gutenberg University) gGmbH. 
                       The project was kindly supported by Intel´s Pandemic Response Technology Initiative.
                       """),
                   html.Br(),
                   html.P("We value your feedback!"),
                   html.P("If you found an error or have a feature request for the dashboard, please, report it:"),
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
                   html.P([
                       "Follow #CoVigator in Twitter ",
                       html.A(html.Img(src="/assets/2021_twitter_logo_blue.png", height="25px"),
                              href="https://twitter.com/hashtag/CoVigator",
                              target="_blank")
                   ])

               ],
               style={"text-align": "left", "font-size": 16}),
        ])
