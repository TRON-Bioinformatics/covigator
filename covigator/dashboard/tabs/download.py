import os
from urllib.parse import quote as urlquote
import dash_html_components as html
from flask import send_from_directory
import dash_bootstrap_components as dbc
import dash_core_components as dcc


def get_tab_download(content_folder):
    return dbc.Card(
        dbc.CardBody(
            children=[
                dcc.Markdown("""
                    ** Download the raw CoVigator data**
                     
                    Every file has a suffix with the version of the form `_v12` and corresponds to a table in our database. 
                    
                    * `variant_observation` contains the variant calls from both ENA and GISAID
                    * `subclonal_variant_observation` contains the variant calls from ENA with a VAF < 80 %
                    * `variant` contains the unique variants without any sample specific information
                    * `variant_cooccurrence` contains the cooccurrence matrix between ENA clonal variants
                    * `sample_ena` contains the ENA samples metadata
                    * `job_ena` contains CoVigator processing metadata on the ENA samples (useful to identify excluded samples)
                    * `sample_gisaid` contains the GISAID samples metadata
                    * `job_gisaid` contains CoVigator processing metadata on the GISAID samples (useful to identify excluded samples)
                    * `conservation` contains the ConsHMM conservation tracks
                    * `gene` contains the gene annotations as provided by Ensembl
                    
                     """, style={"font-size": 16}),
                html.Br(),
                get_downloadable_files(content_folder=content_folder)
            ])
    )


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
                [dbc.ListGroupItem(filename, href="/download/{}".format(urlquote(filename))) for filename in files]))


def set_callbacks_download_tab(app, content_folder):

    @app.server.route("/download/<path:path>")
    def download(path):
        """Serve a file from the upload directory."""
        return send_from_directory(content_folder, path, as_attachment=True)
