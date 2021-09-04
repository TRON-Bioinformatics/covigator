from dask.distributed import Client
from covigator.database.model import Log, DataSource, CovigatorModule
from covigator.processor.ena_processor import EnaProcessor
from covigator.processor.gisaid_processor import GisaidProcessor
from covigator.tests.unit_tests.abstract_test import AbstractTest


class ProcessorTests(AbstractTest):

    def setUp(self) -> None:
        self.ena_processor = EnaProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.gisaid_processor = GisaidProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)

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
