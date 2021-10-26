import dash_html_components as html
import dash_core_components as dcc

import covigator


def get_footer():
    tron_logo = "/assets/tron_logo_no_bg.png"
    return html.Footer(
        [
            html.Div(
                children=[
                    html.Br(),
                    html.P("CoVigator {} © 2021 TRON. All Rights Reserved".format(covigator.VERSION)),
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
                    dcc.Markdown([
                        """
We gratefully acknowledge all data contributors, i.e. the Authors and their Originating laboratories responsible for 
obtaining the specimens, and their Submitting laboratories for generating the genetic sequence and metadata and sharing 
via the GISAID Initiative (1) and the European Nucleotide Archive (2), on which this research is based.

1) Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global 
   health. Global Challenges, 1:33-46. DOI: https://doi.org/10.1002/gch2.1018 PMCID: 31565258

2) Leinonen, R., Akhtar, R., Birney, E., Bower, L., Cerdeno-Tárraga, A., Cheng, Y., Cleland, I., Faruque, N., 
   Goodgame, N., Gibson, R., Hoad, G., Jang, M., Pakseresht, N., Plaister, S., Radhakrishnan, R., Reddy, K., 
   Sobhany, S., Hoopen, P. Ten, Vaughan, R., Zalunin V., Cochrane, G. (2011). The European nucleotide archive. 
   Nucleic Acids Research, 39(SUPPL. 1), D28. https://doi.org/10.1093/nar/gkq967
                        """
                    ]),
                    html.P([
                        html.A("DATA PROTECTION", href="https://tron-mainz.de/data-protection/",  target="_blank"), " | ",
                        html.A("IMPRINT", href="https://tron-mainz.de/imprint/",  target="_blank")])
                ],
                # this bit makes sure the footer sticks at the bottom
                #style={"position": "relative", "bottom": "0", "width": "100%", "overflow": "hidden", "height": "120px"}
            )
        ])
