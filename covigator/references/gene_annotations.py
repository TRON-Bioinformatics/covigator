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
        logger.info(self.config)
        self.session = session
        assert self.config.reference_gene_annotations is not None and \
               os.path.exists(self.config.reference_gene_annotations), \
            "Please configure the gene annotations in the variable {}".format(
                self.config.ENV_COVIGATOR_GENE_ANNOTATIONS)
        assert self.config.reference_proteome is not None and \
               os.path.exists(self.config.reference_proteome), \
            "Please configure the reference proteome in the variable {}".format(
                self.config.ENV_COVIGATOR_REF_PEPTIDE_FASTA)

    def load_data(self):
        # reads the JSON
        with open(self.config.reference_gene_annotations) as fd:
            data = json.load(fd)

        # reads the FASTA
        protein_sequences = {}
        for record in SeqIO.parse(self.config.reference_proteome, "fasta"):
            protein_name = record.description.split("gene_symbol:")[1].split(" ")[0]
            protein_sequences[protein_name] = str(record.seq)

        genomic_sequences = {}
        for record in SeqIO.parse(self.config.reference_genome, "fasta"):
            name = record.description.split(" ")[0]
            genomic_sequences[name] = str(record.seq)

        for seq_name in genomic_sequences:
            self.session.add(Gene(identifier=seq_name, name=seq_name, sequence=genomic_sequences[seq_name]))

        for g in data["genes"]:
            self.session.add(Gene(identifier=g["id"], name=g["name"], start=int(g["start"]), end=int(g["end"]),
                                  sequence=protein_sequences.get(g["name"]), data=g))
        self.session.commit()
        logger.info("Loaded into the database {} genes".format(len(data["genes"])))
