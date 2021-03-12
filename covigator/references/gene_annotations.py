import shutil
from contextlib import closing
from urllib import request
import os
import json
from sqlalchemy.orm import Session
from covigator import ENV_COVIGATOR_STORAGE_FOLDER
from covigator.model import Gene


class GeneAnnotationsLoader:

    def __init__(self, session: Session):
        self.storage_folder = os.getenv(ENV_COVIGATOR_STORAGE_FOLDER, "./data/covigator")
        # TODO: this is hardcoded which fixes covigator to Sars-Cov-2, but the problem is that other infectious
        #  organisms do not have their reference in JSON format, this is new and more complete than GFF
        self.reference_file = "ftp://ftp.ensemblgenomes.org/pub/viruses/json/sars_cov_2/sars_cov_2.json"
        self.session = session

    def load_data(self):
        file_name = os.path.join(self.storage_folder, self.reference_file.split('/')[-1])
        # download JSON
        with closing(request.urlopen(self.reference_file)) as r:
            with open(file_name, 'wb') as f:
                shutil.copyfileobj(r, f)
        # reads the JSON
        data = json.load(open(file_name))
        for g in data["genes"]:
            self.session.add(Gene(identifier=g["id"], name=g["name"], data=g))
        self.session.commit()
