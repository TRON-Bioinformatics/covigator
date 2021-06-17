import os
import pkg_resources
import covigator
from covigator.configuration import Configuration


class FakeConfiguration(Configuration):

    def __init__(self):
        # use a folder with write permissions for testing
        os.environ[self.ENV_COVIGATOR_STORAGE_FOLDER] = "./data/covigator"
        os.environ[self.ENV_COVIGATOR_REF_FASTA] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/MN908947.3.fa")
        os.environ[self.ENV_COVIGATOR_GENE_ANNOTATIONS] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/sars_cov_2.json")
        super().__init__()
