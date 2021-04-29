from unittest import TestCase
from dask.distributed import Client
from covigator.database.database import Database
from covigator.database.model import Log, DataSource, CovigatorModule
from covigator.processor.ena_processor import EnaProcessor
from faker import Faker
from covigator.processor.gisaid_processor import GisaidProcessor


class ProcessorTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()
        self.ena_processor = EnaProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1))
        self.gisaid_processor = GisaidProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1))
        self.faker = Faker()

    def test_no_ena_jobs(self):
        self.ena_processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.ENA)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)

    def test_no_gisaid_jobs(self):
        self.gisaid_processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.GISAID)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)
