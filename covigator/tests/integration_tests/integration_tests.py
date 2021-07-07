from unittest import TestCase
from dask.distributed import Client
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.processor.ena_processor import EnaProcessor
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID
from covigator.dashboard.dashboard import Dashboard
import time


class IntegrationTests(TestCase):

    def setUp(self) -> None:
        self.config = Configuration()

    def test_access(self):
        """
        If given enough time, this creates a database with all runs complying with the selection criteria from ENA
        """
        accessor = EnaAccessor(tax_id=SARS_COV_2_TAXID, host_tax_id=HOMO_SAPIENS_TAXID,
                               database=Database(config=self.config))
        accessor.access()

    def test_process_dask(self):
        """
        If given enough time, this processes all ENA samples in the database
        """
        client = Client(n_workers=6, threads_per_worker=1)
        processor = EnaProcessor(database=Database(), dask_client=client, config=self.config)
        processor.process()

    def test_dashboard(self):
        dashboard = Dashboard(config=self.config)
        dashboard.start_dashboard(debug=True)
        time.sleep(60)
