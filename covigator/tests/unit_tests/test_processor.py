from unittest import TestCase

from dask.distributed import Client

from covigator.database.database import Database
from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna, JobStatus, JobEna
from covigator.processor.dummy_processor import DummyProcessor
from covigator.processor.ena_processor import EnaProcessor
from faker import Faker

from covigator.processor.gisaid_processor import GisaidProcessor
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample


class ProcessorTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()
        self.ena_processor = EnaProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), batch_size=10)
        self.gisaid_processor = GisaidProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), batch_size=10)
        self.dummy_processor = DummyProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), batch_size=10)
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
        self.assertEqual(data.get("batches"), 0)

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
        self.assertEqual(data.get("batches"), 0)

    def test_dummy_processor_no_jobs(self):
        self.dummy_processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.ENA)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)
        self.assertEqual(data.get("batches"), 0)

    def test_dummy_processor(self):

        for _ in range(50):
            sample_ena, sample, job = get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.PENDING)
            self.session.add(sample_ena)
            self.session.commit()
            self.session.add_all([sample, job])
            self.session.commit()

        self.dummy_processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.ENA)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 50)
        self.assertFalse(log.has_error)
        data = log.data
        self.assertEqual(data.get("processed"), 50)
        self.assertEqual(data.get("batches"), 5)
        self.assertEqual(self.session.query(JobEna).filter(JobEna.status == JobStatus.COOCCURRENCE).count(), 50)
