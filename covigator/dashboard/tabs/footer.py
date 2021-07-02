import dash_html_components as html

import covigator


def get_footer():
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
