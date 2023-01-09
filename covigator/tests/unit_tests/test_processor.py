import covigator
import pkg_resources
from dask.distributed import Client
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus, SampleEna
from covigator.processor.ena_processor import EnaProcessor
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.faked_objects import FakeEnaProcessor, \
    FakeProcessorFailing
from covigator.tests.unit_tests.mocked import mock_samples


class ProcessorTests(AbstractTest):

    def setUp(self) -> None:

        dask_client = Client(n_workers=int(1), threads_per_worker=1)

        self.processor = EnaProcessor(
            database=self.database, dask_client=dask_client, config=self.config, wait_time=1)
        self.fake_processor = FakeEnaProcessor(
            database=self.database, dask_client=dask_client, config=self.config)
        self.fake_processor_fail = FakeProcessorFailing(
            database=self.database, dask_client=dask_client, config=self.config, source=DataSource.ENA)

    def test_no_jobs(self):
        self.processor.process()
        self.assertEqual(self.session.query(Log).count(), 1)
        log = self.session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.module, CovigatorModule.PROCESSOR)
        self.assertEqual(log.processed, 0)
        data = log.data
        self.assertEqual(data.get("processed"), 0)

    def test_fake_processor(self):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.DOWNLOADED)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.DOWNLOADED), 10)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.FINISHED), 0)
        self.fake_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=DataSource.ENA, status=JobStatus.FINISHED), 10)
        finished_jobs = self.queries.find_first_by_status(
            data_source=DataSource.ENA, status=(JobStatus.FINISHED, ), n=10)
        for j in finished_jobs:
            self.assertEqual(j.status, JobStatus.FINISHED)
            self.assertIsNotNone(j.analysed_at)
            self.assertIsNotNone(j.pangolin_lineage)

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
        self.fake_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 1)

    def test_does_not_process_already_queued(self, source=DataSource.ENA):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.QUEUED,
                     source=source.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.QUEUED), 10)
        self.fake_processor.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.PENDING), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.QUEUED), 10)

    def test_failed_sample(self, source=DataSource.ENA):
        mock_samples(faker=self.faker, session=self.session, num_samples=10, job_status=JobStatus.DOWNLOADED,
                     source=source.name)

        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 10)
        self.fake_processor_fail.process()
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.DOWNLOADED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FINISHED), 0)
        self.assertEqual(self.queries.count_jobs_by_status(data_source=source, status=JobStatus.FAILED_PROCESSING), 10)

    def test_load_pangolin_with_none(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.with_none.csv")
        sample = SampleEna(run_accession="TEST", fastq_ftp="", fastq_md5="", num_fastqs=1)
        EnaProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_lineage, "")

        self.session.add(sample)
        self.session.commit()

    def test_load_pangolin(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.pangolin.csv")
        sample = SampleEna(run_accession="TEST", fastq_ftp="blabla", fastq_md5="blabla", num_fastqs=1)
        sample = EnaProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_pangolin_version, "4.1.2")
        self.assertEqual(sample.pangolin_version, "PANGO-v1.14")
        self.assertEqual(sample.pangolin_scorpio_version, "0.3.17")
        self.assertEqual(sample.pangolin_constellation_version, "v0.1.10")
        self.assertEqual(sample.pangolin_qc_status, "pass")
        self.assertEqual(sample.pangolin_note, "Assigned from designation hash.")
        self.assertEqual(sample.pangolin_lineage, "B")
        self.assertEqual(sample.pangolin_conflict, 0.0)
        self.assertEqual(sample.pangolin_ambiguity_score, 0.0)
        self.session.add(sample)
        self.session.commit()

    def test_load_another_pangolin(self):
        pangolin_path = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.another_pangolin.csv")
        sample = SampleEna(run_accession="TEST", fastq_ftp="blabla", fastq_md5="blabla", num_fastqs=1)
        sample = EnaProcessor.load_pangolin(sample, path=pangolin_path)
        self.assertEqual(sample.pangolin_pangolin_version, "4.1.2")
        self.assertEqual(sample.pangolin_version, "PUSHER-v1.14")
        self.assertEqual(sample.pangolin_scorpio_version, "0.3.17")
        self.assertEqual(sample.pangolin_constellation_version, "v0.1.10")
        self.assertEqual(sample.pangolin_qc_status, "pass")
        self.assertEqual(sample.pangolin_note, "Usher placements: AY.44(2/2)")
        self.assertEqual(sample.pangolin_lineage, "AY.44")
        self.assertEqual(sample.pangolin_conflict, 0.0)
        self.assertEqual(sample.pangolin_ambiguity_score, 0.0)
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

