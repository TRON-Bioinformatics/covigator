from unittest import TestCase

from dask.distributed import Client

from covigator.database.database import Database
from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna
from covigator.processor.ena_processor import EnaProcessor
from faker import Faker


class EnaProcessorTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()
        self.processor = EnaProcessor(database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1))
        self.faker = Faker()

    def test_no_jobs(self):
        self.processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.ENA)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)
