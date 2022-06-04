from dash import html
from dash import dcc
import dash_bootstrap_components as dbc

import covigator


def get_footer():
    tron_logo = "/assets/tron_logo_no_bg.png"
    return html.Footer(
        [
            dbc.Row([
                dbc.Col([
                    html.Br(),
                    html.P("CoVigator {} © 2021 TRON. All Rights Reserved".format(covigator.VERSION)),
                    html.P([
                        "GISAID data provided on this website are subject to GISAID’s ",
                        html.A("Terms and Conditions", href="https://www.gisaid.org/registration/terms-of-use/",
                               target="_blank")]),
                    html.A(html.Img(src=tron_logo, id="tron-logo"), href="https://tron-mainz.de",  target="_blank"),
                    html.Br(),
                    dcc.Markdown("""
                    TRON is an independent biopharmaceutical non-profit translational research organization pursuing 
                    new diagnostics and drugs for the treatment of cancer and other diseases with high medical need. 
                    We focus our transdisciplinary competencies in genomics and immunology to 
                    1) develop novel platforms for the identification and validation of “omics”-based biomarkers
                    and 2) for harnessing and modulating immune system components, for use in personalized therapies.
                    Partnering with academia and industry, TRON executes research at the leading edge to support 
                    innovative drug design for human health.
                    """),
                    html.P([
                        """
                    Intel is committed to accelerating access to technology that can combat the current pandemic 
                    and enable scientific discovery that better prepares our world for future crises. Funding for this 
                    solution was funded in part by """,
                        html.A(
                            "Intel’s Pandemic Response Technology Initiative",
                            href="https://newsroom.intel.com/news/intel-commits-technology-response-combat-coronavirus/",
                            target="_blank"),
                        """. For more 
                    information about healthcare solutions from Intel, visit intel.com/healthcare. For more 
                    information about Intel’s COVID-19 response, visit """,
                        html.A(
                            "intel.com/COVID-19",
                            href="https://www.intel.com/content/www/us/en/corporate-responsibility/covid-19-response.html",
                            target="_blank"),
                        "."
                        ]),
                    html.P([
                        html.A("DATA PROTECTION", href="https://tron-mainz.de/data-protection/",  target="_blank"), " | ",
                        html.A("IMPRINT", href="https://tron-mainz.de/imprint/",  target="_blank")])
                ], style={"color": "#003c78", 'margin-left': '15px'})]
                # this bit makes sure the footer sticks at the bottom
                #style={"position": "relative", "bottom": "0", "width": "100%", "overflow": "hidden", "height": "120px"}
            )
        ])
