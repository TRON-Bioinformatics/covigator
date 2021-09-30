import os
import pkg_resources
import covigator
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID


class FakeConfiguration(Configuration):

    def __init__(self):
        # use a folder with write permissions for testing
        os.environ[self.ENV_COVIGATOR_STORAGE_FOLDER] = "./data/covigator"
        os.environ[self.ENV_COVIGATOR_REF_FASTA] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/MN908947.3.fa")
        os.environ[self.ENV_COVIGATOR_GENE_ANNOTATIONS] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/sars_cov_2.json")
        os.environ[self.ENV_COVIGATOR_GENE_DN_DS_ANNOTATIONS] = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/sars_cov_2_dn_ds.json")

        # this makes sure that we do not wipe a relevant database by mistake
        os.environ[self.ENV_COVIGATOR_TABLE_VERSION] = "_test"

        super().__init__()


class FakeEnaAccessor(EnaAccessor):

    def __init__(self, results, database=None):
        # uses an in memory database or the one provided
        super().__init__(tax_id=SARS_COV_2_TAXID, host_tax_id=HOMO_SAPIENS_TAXID,
                         database=database if database else Database(test=True, config=Configuration()))
        self.results = results

    def _get_ena_runs_page(self):
        return self.results
