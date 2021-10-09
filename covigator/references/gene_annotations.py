import os
import json

import pandas as pd
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
        assert self.config.reference_gene_ns_s_fractions is not None and os.path.exists(
            self.config.reference_gene_ns_s_fractions), \
            "Please configure the gene annotations for NS and S fractions in the variable {}".format(
                self.config.ENV_COVIGATOR_GENE_NS_S_FRACTIONS)
        assert self.config.reference_domain_ns_s_fractions is not None and os.path.exists(
            self.config.reference_domain_ns_s_fractions), \
            "Please configure the domain annotations for NS and S fractions in the variable {}".format(
                self.config.ENV_COVIGATOR_DOMAIN_NS_S_FRACTIONS)

    def load_data(self):
        # reads the JSON
        with open(self.config.reference_gene_annotations) as fd:
            data = json.load(fd)

        # NOTE: load fraction of synonymous and non synonymous
        genes_fractions = pd.read_csv(self.config.reference_gene_ns_s_fractions)
        domains_fractions = pd.read_csv(self.config.reference_domain_ns_s_fractions)

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
