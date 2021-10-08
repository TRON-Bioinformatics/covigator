from dask.distributed import Client
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus, JobEna
from covigator.processor.ena_processor import EnaProcessor
from covigator.processor.gisaid_processor import GisaidProcessor
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.faked_objects import FakeEnaProcessor, FakeEnaProcessorExcludingSamples, \
    FakeGisaidProcessor, FakeEnaProcessorFailing
from covigator.tests.unit_tests.mocked import mock_samples


class ProcessorTests(AbstractTest):

    def setUp(self) -> None:

        self.ena_processor = EnaProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.gisaid_processor = GisaidProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.fake_ena_processor = FakeEnaProcessor(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.fake_ena_processor_excluder = FakeEnaProcessorExcludingSamples(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.fake_ena_processor_fail = FakeEnaProcessorFailing(
            database=self.database, dask_client=Client(n_workers=int(1), threads_per_worker=1), config=self.config)
        self.fake_gisaid_processor = FakeGisaidProcessor(
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

    def test_ena_processor(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.PENDING,
                     source=DataSource.ENA.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 10)
        self.fake_ena_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.FINISHED), 10)

    def test_excluded_sample(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.PENDING,
                     source=DataSource.ENA.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 10)
        self.fake_ena_processor_excluder.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.FINISHED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.EXCLUDED), 10)

    def test_gisaid_processor(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.PENDING,
                     source=DataSource.GISAID.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.GISAID, status=JobStatus.PENDING), 10)
        self.fake_gisaid_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.GISAID, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.GISAID, status=JobStatus.FINISHED), 10)

    def test_does_not_process_already_queued(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.QUEUED,
                     source=DataSource.ENA.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.QUEUED), 10)
        self.fake_ena_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.QUEUED), 10)

    def test_failed_sample(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.PENDING,
                     source=DataSource.ENA.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 10)
        self.fake_ena_processor_fail.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.FINISHED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA,
                                                           status=JobStatus.FAILED_PROCESSING), 10)
