import os
import json
from Bio import SeqIO
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.database.model import Gene, Domain
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

        # NOTE: load fraction of synonymous and non synonymous
        with open(self.config.reference_gene_dn_ds_annotations) as fd:
            data_dn_ds = json.load(fd)

        count_genes = 0
        count_domains = 0

        for g in data["genes"]:
            # persist a gene
            gene = Gene(identifier=g["id"], name=g["name"], start=int(g["start"]), end=int(g["end"]),
                        fraction_synonymous=0.0, fraction_non_synonymous=0.0)
            self.session.add(gene)
            self.session.commit()
            count_genes += 1

            for d in g.get("transcripts", [])[0].get("translations", [])[0].get("protein_features"):
                if d.get("dbname") == "Pfam":
                    try:
                        domain = Domain(name=d["description"], start=int(d["start"]), end=int(d["end"]),
                                        fraction_synonymous=0.0, fraction_non_synonymous=0.0,
                                        gene_identifier=gene.identifier, gene_name=gene.name)
                        self.session.add(domain)
                        self.session.commit()
                    except IntegrityError:
                        logger.warn("Domain {} is duplicated in input with coordinates [{}, {}]".format(
                            d["description"], d["start"], d["end"]))
                        self.session.rollback()
                    count_domains += 1

        logger.info("Loaded into the database {} genes and {} domains".format(count_genes, count_domains))

    def get_dn_ds_by_gene(self, data, gene):
        dnds = None
        for e in data:
            if e.get("name") == gene:
                dnds = e.get("ratio_synonymous_non_synonymous")
                break
        return dnds
