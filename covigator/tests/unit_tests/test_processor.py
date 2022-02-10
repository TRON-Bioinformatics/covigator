import numpy as np

import covigator
import pkg_resources
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

    def test_load_pangolin(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.csv")
        job = JobEna(lofreq_pangolin_path=pangolin_path)
        EnaProcessor.load_lofreq_pangolin(job)
        self.assertEqual(job.pangolin_pangolin_version, "3.1.19")
        self.assertEqual(job.pangolin_version, "PLEARN-v1.2.123")
        self.assertEqual(job.pangolin_pangoLEARN_version, "2022-01-20")
        self.assertEqual(job.pangolin_pango_version, "v1.2.123")
        self.assertEqual(job.pangolin_status, "passed_qc")
        self.assertEqual(job.pangolin_note, "")
        self.assertEqual(job.pangolin_lineage, "B.1.336")
        self.assertEqual(job.pangolin_conflict, 0.0)
        self.assertEqual(job.pangolin_ambiguity_score, 1.0)

    def test_load_dedup_metrics(self):
        deduplication_metrics_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.deduplication_metrics.txt")
        job = JobEna(deduplication_metrics_path=deduplication_metrics_path)
        EnaProcessor.load_deduplication_metrics(job)
        self.assertEqual(job.percent_duplication, 0.207048)
        self.assertEqual(job.unpaired_reads_examined, 2497)
        self.assertEqual(job.read_pairs_examined, 0)
        self.assertEqual(job.secondary_or_supplementary_reads, 7)
        self.assertEqual(job.unmapped_reads, 3)
        self.assertEqual(job.unpaired_read_duplicates, 517)
        self.assertEqual(job.read_pair_duplicates, 0)
        self.assertEqual(job.read_pair_optical_duplicates, 0)

    def test_load_horizontal_coverage(self):
        horizontal_coverage_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.coverage.tsv")
        job = JobEna(horizontal_coverage_path=horizontal_coverage_path)
        EnaProcessor.load_coverage_results(job)
        self.assertEqual(job.mean_depth, 8.19399)
        self.assertEqual(job.mean_base_quality, 35.9)
        self.assertEqual(job.mean_mapping_quality, 60.0)
        self.assertEqual(job.num_reads, 1984)
        self.assertEqual(job.covered_bases, 26057)
        self.assertEqual(job.coverage, 87.1384)

