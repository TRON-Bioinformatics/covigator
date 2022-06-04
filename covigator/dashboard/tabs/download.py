import os
from urllib.parse import quote as urlquote
from dash import html
from flask import send_from_directory
import dash_bootstrap_components as dbc
from dash import dcc
from logzero import logger


def get_tab_download(content_folder):
    return dbc.CardBody(
            children=[
                dcc.Markdown("""
                    ** Download the raw CoVigator data derived from ENA** 
                    
                    * `variant_observation` contains the variant calls from ENA
                    * `subclonal_variant_observation` contains the variant calls from ENA with a VAF < 80 % and >= 5 %
                    * `low_frequency_variant_observation` contains the variant calls from ENA with a VAF < 5 %
                    * `variant` contains the unique variants without any sample specific information
                    * `variant_cooccurrence` contains the cooccurrence matrix between ENA clonal variants
                    * `sample_ena` contains the ENA samples metadata and some sample level derived data (eg: pangolin, coverage analysis, etc.)
                    * `conservation` contains the ConsHMM conservation tracks
                    * `gene` contains the gene annotations as provided by Ensembl
                    * `domain` contains the Pfam protein domains
                    
                    **NOTE**: no GISAID data, metadata or derived data that enables reverse engineer GISAID original 
                    sequences is available for download in agreement with GISAID's terms and conditions 
                     """, style={"font-size": 16}),
                html.Br(),
                get_downloadable_files(content_folder=content_folder)
            ])


def get_downloadable_files(content_folder):
    files = []
    if content_folder:
        for filename in os.listdir(content_folder):
            path = os.path.join(content_folder, filename)
            if os.path.isfile(path):
                files.append(filename)
    if len(files) == 0:
        return [html.Li("No files to download yet!")]
    else:
        return html.Div(
            children=dbc.ListGroup(
                [dbc.ListGroupItem(html.A(filename, href="/download/{}".format(urlquote(filename)),  target="_blank"))
                 for filename in files]))


def set_callbacks_download_tab(app, content_folder):

    @app.server.route("/download/<path:path>")
    def download(path):
        """Serve a file from the upload directory."""
        logger.info("hey")
        return send_from_directory(content_folder, path, as_attachment=True)
