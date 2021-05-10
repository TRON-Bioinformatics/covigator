import os

from covigator.configuration import Configuration


class FakeConfiguration(Configuration):

    def __init__(self):
        # use a folder with write permissions for testing
        os.environ[self.ENV_COVIGATOR_STORAGE_FOLDER] = "./data/covigator"
        os.environ[self.ENV_COVIGATOR_REF_FASTA] = "../resources/MN908947.3.fa"
        os.environ[self.ENV_COVIGATOR_GENE_ANNOTATIONS] = "../resources/sars_cov_2.json"
        super().__init__()
