import os

from covigator.configuration import Configuration


class FakeConfiguration(Configuration):

    def __init__(self):
        # use a folder with write permissions for testing
        os.environ[self.ENV_COVIGATOR_STORAGE_FOLDER] = "./data/covigator"
        os.environ[self.ENV_COVIGATOR_REF_FASTA] = "../resources/MN908947.3.fa"
        os.environ[self.ENV_COVIGATOR_REF_PEPTIDE_FASTA] = "../resources/Sars_cov_2.ASM985889v3.pep.all.fa"
        os.environ[self.ENV_COVIGATOR_GENE_ANNOTATIONS] = "../resources/sars_cov_2.json"
        super().__init__()
