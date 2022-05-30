import numpy as np
from parameterized import parameterized

import covigator
import pkg_resources
from dask.distributed import Client
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus, SampleEna, SampleGisaid
from covigator.processor.ena_processor import EnaProcessor
from covigator.processor.gisaid_processor import GisaidProcessor
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.faked_objects import FakeEnaProcessor, \
    FakeGisaidProcessor, FakeProcessorFailing
from covigator.tests.unit_tests.mocked import mock_samples


class ProcessorTests(AbstractTest):

    def setUp(self) -> None:

        dask_client = Client(n_workers=int(1), threads_per_worker=1)

        ena_processor = EnaProcessor(
            database=self.database, dask_client=dask_client, config=self.config, wait_time=1)
        gisaid_processor = GisaidProcessor(
            database=self.database, dask_client=dask_client, config=self.config, wait_time=1)

        self.processors = {
            DataSource.ENA: ena_processor,
            DataSource.GISAID: gisaid_processor
        }

        fake_ena_processor = FakeEnaProcessor(
            database=self.database, dask_client=dask_client, config=self.config)
        fake_gisaid_processor = FakeGisaidProcessor(
            database=self.database, dask_client=dask_client, config=self.config)
        self.fake_processors = {
            DataSource.ENA: fake_ena_processor,
            DataSource.GISAID: fake_gisaid_processor
        }

        fake_ena_processor_fail = FakeProcessorFailing(
            database=self.database, dask_client=dask_client, config=self.config, source=DataSource.ENA)
        fake_gisaid_processor_fail = FakeProcessorFailing(
            database=self.database, dask_client=dask_client, config=self.config, source=DataSource.GISAID)
        self.failing_processors = {
            DataSource.ENA: fake_ena_processor_fail,
            DataSource.GISAID: fake_gisaid_processor_fail
        }

    @parameterized.expand([(DataSource.ENA, ), (DataSource.GISAID, )])
    def test_no_jobs(self, source):
        self.processors.get(source).process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, source)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)

    @parameterized.expand([(DataSource.ENA, ), (DataSource.GISAID, )])
    def test_fake_processor(self, source):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.DOWNLOADED,
                     source=source.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 10)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 0)
        self.fake_processors.get(source).process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 10)
        finished_jobs = self.queries.find_first_by_status(data_source=source, status=[JobStatus.FINISHED], n=10)
        for j in finished_jobs:
            self.assertEqual(j.status, JobStatus.FINISHED)
            self.assertIsNotNone(j.analysed_at)
            self.assertIsNotNone(j.pangolin_lineage)

    # @parameterized.expand([(DataSource.ENA,), (DataSource.GISAID,)])
    def test_processor(self, source=DataSource.ENA):
        samples = mock_samples(faker=self.faker, session=self.session, num_samples=1, job_status=JobStatus.DOWNLOADED,
                     source=source.name)
        sample = samples[0]

        fastq1 = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test_data_1.fastq.gz")
        fastq2 = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test_data_2.fastq.gz")
        sample.fastq_path = "{fastq1},{fastq2}".format(fastq1=fastq1, fastq2=fastq2)
        self.session.merge(sample)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 1)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 0)
        self.fake_processors.get(source).process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 1)

    @parameterized.expand([(DataSource.ENA, ), (DataSource.GISAID, )])
    def test_does_not_process_already_queued(self, source):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.QUEUED,
                     source=source.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.QUEUED), 10)
        self.fake_processors.get(source).process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.QUEUED), 10)

    @parameterized.expand([(DataSource.ENA, ), (DataSource.GISAID, )])
    def test_failed_sample(self, source):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.DOWNLOADED,
                     source=source.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 10)
        self.failing_processors.get(source).process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FAILED_PROCESSING), 10)

    def test_load_pangolin_from_ena(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.csv")
        sample = SampleEna(run_accession="TEST", fastq_ftp="", fastq_md5="", num_fastqs=1)
        EnaProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_pangolin_version, "3.1.19")
        self.assertEqual(sample.pangolin_version, "PLEARN-v1.2.123")
        self.assertEqual(sample.pangolin_pangoLEARN_version, "2022-01-20")
        self.assertEqual(sample.pangolin_pango_version, "v1.2.123")
        self.assertEqual(sample.pangolin_status, "passed_qc")
        self.assertEqual(sample.pangolin_note, "")
        self.assertEqual(sample.pangolin_lineage, "B.1.336")
        self.assertEqual(sample.pangolin_conflict, 0.0)
        self.assertEqual(sample.pangolin_ambiguity_score, 1.0)

        self.session.add(sample)
        self.session.commit()

    def test_load_pangolin_with_none(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.with_none.csv")
        sample = SampleEna(run_accession="TEST", fastq_ftp="", fastq_md5="", num_fastqs=1)
        EnaProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_lineage, "")

        self.session.add(sample)
        self.session.commit()

    def test_load_pangolin_from_gisaid(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.csv")
        sample = SampleGisaid(run_accession="TEST")
        GisaidProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_pangolin_version, "3.1.19")
        self.assertEqual(sample.pangolin_version, "PLEARN-v1.2.123")
        self.assertEqual(sample.pangolin_pangoLEARN_version, "2022-01-20")
        self.assertEqual(sample.pangolin_pango_version, "v1.2.123")
        self.assertEqual(sample.pangolin_status, "passed_qc")
        self.assertEqual(sample.pangolin_note, "")
        self.assertEqual(sample.pangolin_lineage, "B.1.336")
        self.assertEqual(sample.pangolin_conflict, 0.0)
        self.assertEqual(sample.pangolin_ambiguity_score, 1.0)

        self.session.add(sample)
        self.session.commit()

    def test_load_dedup_metrics(self):
        deduplication_metrics_path = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/test.deduplication_metrics.txt")
        sample = SampleEna(run_accession="TEST", deduplication_metrics_path=deduplication_metrics_path, fastq_ftp="", fastq_md5="",
                           num_fastqs=1)
        EnaProcessor.load_deduplication_metrics(sample)
        self.assertEqual(sample.percent_duplication, 0.919339)
        self.assertEqual(sample.unpaired_reads_examined, 0)
        self.assertEqual(sample.read_pairs_examined, 352625)
        self.assertEqual(sample.secondary_or_supplementary_reads, 1972)
        self.assertEqual(sample.unmapped_reads, 0)
        self.assertEqual(sample.unpaired_read_duplicates, 0)
        self.assertEqual(sample.read_pair_duplicates, 324182)
        self.assertEqual(sample.read_pair_optical_duplicates, 0)

        self.session.add(sample)
        self.session.commit()

    def test_load_horizontal_coverage(self):
        horizontal_coverage_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.coverage.tsv")
        sample = SampleEna(run_accession="TEST", horizontal_coverage_path=horizontal_coverage_path, fastq_ftp="",
                           fastq_md5="", num_fastqs=1)
        EnaProcessor.load_coverage_results(sample)
        self.assertEqual(sample.mean_depth, 8.19399)
        self.assertEqual(sample.mean_base_quality, 35.9)
        self.assertEqual(sample.mean_mapping_quality, 60.0)
        self.assertEqual(sample.num_reads, 1984)
        self.assertEqual(sample.covered_bases, 26057)
        self.assertEqual(sample.coverage, 87.1384)

        self.session.add(sample)
        self.session.commit()

