from unittest import TestCase
from dask.distributed import Client
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.database import Database
from covigator.processor.processor import Processor
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID
import covigator.dashboard.dashboard as dashboard
import time


class IntegrationTests(TestCase):

    def test_access(self):
        """
        If given enough time, this creates a database with all runs complying with the selection criteria from ENA
        """
        accessor = EnaAccessor(tax_id=SARS_COV_2_TAXID, host_tax_id=HOMO_SAPIENS_TAXID, database=Database())
        accessor.access()

    def test_process_dask(self):

        client = Client(n_workers=6, threads_per_worker=1)
        processor = Processor(database=Database(), dask_client=client)
        processor.process()

    def test_dashboard(self):
        dashboard.main()
        time.sleep(60)
