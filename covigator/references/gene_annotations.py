import os
import json
from Bio import SeqIO
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.database.model import Gene
from logzero import logger


class GeneAnnotationsLoader:

    def __init__(self, session: Session, config: Configuration):
        self.config = config
        self.session = session
        assert self.config.reference_gene_annotations is not None and os.path.exists(
            self.config.reference_gene_annotations), \
            "Please configure the gene annotations in the variable {}".format(self.config.ENV_COVIGATOR_GENE_ANNOTATIONS)
        assert self.config.reference_gene_dn_ds_annotations is not None and os.path.exists(
            self.config.reference_gene_dn_ds_annotations), \
            "Please configure the gene annotations for dN/dS in the variable {}".format(
                self.config.ENV_COVIGATOR_GENE_DN_DS_ANNOTATIONS)

    def load_data(self):
        # reads the JSON
        with open(self.config.reference_gene_annotations) as fd:
            data = json.load(fd)

        with open(self.config.reference_gene_dn_ds_annotations) as fd:
            data_dn_ds = json.load(fd)

        for g in data["genes"]:
            self.session.add(Gene(identifier=g["id"], name=g["name"], start=int(g["start"]), end=int(g["end"]),
                                  ratio_synonymous_non_synonymous=self.get_dn_ds_by_gene(data_dn_ds, g["name"]),
                                  data=g))

        self.session.commit()
        logger.info("Loaded into the database {} genes".format(len(data["genes"])))

    def get_dn_ds_by_gene(self, data, gene):
        dnds = None
        for e in data:
            if e.get("name") == gene:
                dnds = e.get("ratio_synonymous_non_synonymous")
                break
        return dnds
