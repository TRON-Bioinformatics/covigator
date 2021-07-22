import dash_html_components as html

import covigator


def get_footer():
    tron_logo = "/assets/tron_logo_no_bg.png"
    return html.Footer(
        [
            html.Div(
                children=[
                    html.Br(),
                    html.P("CoVigator {} © 2021 TRON. All Rights Reserved".format(covigator.VERSION)),
                    html.A(html.Img(src=tron_logo, id="tron-logo"), href="https://tron-mainz.de"),
                    html.P("The CoVigator project is enabled by Intel Corporation servers."),
                    html.P([
                        html.A("DATA PROTECTION", href="https://tron-mainz.de/data-protection/"), " | ",
                        html.A("IMPRINT", href="https://tron-mainz.de/imprint/")])
                ],
                # this bit makes sure the footer sticks at the bottom
                #style={"position": "relative", "bottom": "0", "width": "100%", "overflow": "hidden", "height": "120px"}
            )
        ])
