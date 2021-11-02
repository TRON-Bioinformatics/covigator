import os
import json

import pandas as pd
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session
from covigator.database.model import Gene, Domain
from logzero import logger


class GeneAnnotationsLoader:

    GENES_NS_S_COUNTS_FILENAME = "genes_NS_S.csv"
    DOMAINS_NS_S_COUNTS_FILENAME = "domains_NS_S.csv"
    GENE_ANNOTATIONS_FILENAME = "sars_cov_2.json"

    def __init__(self, session: Session):
        self.session = session

        self.genes_ns_s_counts = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.GENES_NS_S_COUNTS_FILENAME)
        self.domains_ns_s_counts = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.DOMAINS_NS_S_COUNTS_FILENAME)
        self.gene_annotations = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.GENE_ANNOTATIONS_FILENAME)

    def _remove_duplicated_genes(self, data):
        results = {}
        for g in data["genes"]:
            if g["name"] not in results:
                results[g["name"]] = g
            else:
                if int(g["end"]) - int(g["start"]) > int(results[g["name"]]["end"]) - int(results[g["name"]]["start"]):
                    # in case of duplications it keeps only the longer gene
                    results[g["name"]] = g
        data["genes"] = list(results.values())
        return data

    def load_data(self):
        # reads the JSON
        with open(self.gene_annotations) as fd:
            data = json.load(fd)

        data = self._remove_duplicated_genes(data)

        # NOTE: load fraction of synonymous and non synonymous
        genes_fractions = pd.read_csv(self.genes_ns_s_counts)
        domains_fractions = pd.read_csv(self.domains_ns_s_counts)

        count_genes = 0
        count_domains = 0

        for g in data["genes"]:
            # persist a gene
            transcript_id = g["transcripts"][0]["id"]
            gene = Gene(
                identifier=g["id"],
                name=g["name"],
                start=int(g["start"]),
                end=int(g["end"]),
                fraction_synonymous=genes_fractions[genes_fractions.transcript == transcript_id].S.iloc[0],
                fraction_non_synonymous=genes_fractions[genes_fractions.transcript == transcript_id].NS.iloc[0])
            self.session.add(gene)
            self.session.commit()
            count_genes += 1

            for d in g.get("transcripts", [])[0].get("translations", [])[0].get("protein_features"):
                if d.get("dbname") == "Pfam":
                    try:
                        domain = Domain(
                            name=d["description"],
                            description=d["interpro_description"],
                            start=int(d["start"]),
                            end=int(d["end"]),
                            fraction_synonymous=domains_fractions[domains_fractions.domain == d["description"]].S.iloc[0],
                            fraction_non_synonymous=domains_fractions[domains_fractions.domain == d["description"]].NS.iloc[0],
                            gene_identifier=gene.identifier,
                            gene_name=gene.name)
                        self.session.add(domain)
                        self.session.commit()
                    except IntegrityError:
                        logger.warn("Domain {} is duplicated in input with coordinates [{}, {}]".format(
                            d["description"], d["start"], d["end"]))
                        self.session.rollback()
                    count_domains += 1

        logger.info("Loaded into the database {} genes and {} domains".format(count_genes, count_domains))
