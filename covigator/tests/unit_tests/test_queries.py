from unittest import TestCase, skip
from faker import Faker
from covigator.database.database import Database
from covigator.database.model import JobStatus, DataSource
from covigator.database.queries import get_date_of_first_ena_sample, get_date_of_most_recent_ena_sample, \
    get_date_of_last_check, get_date_of_last_update
from covigator.tests.unit_tests.mocked import get_mocked_ena_sample, get_mocked_log


class QueriesTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()
        self.faker = Faker()

    def test_get_date_of_first_ena_sample(self):
        first_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = [get_mocked_ena_sample(faker=self.faker) for _ in range(50)] + \
                  [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.FAILED_LOAD) for _ in range(50)]
        for sample_ena, sample, job in samples:
            self.session.add_all([sample_ena, sample, job])
            if job.status == JobStatus.LOADED:
                if first_sample_date is None:
                    first_sample_date = sample_ena.first_created
                if sample_ena.first_created < first_sample_date:
                    first_sample_date = sample_ena.first_created
        self.session.commit()
        observed_date = get_date_of_first_ena_sample(session=self.session)
        self.assertEqual(observed_date, first_sample_date.date())

    def test_get_date_of_first_ena_sample_empty(self):
        observed_date = get_date_of_first_ena_sample(session=self.session)
        self.assertIsNone(observed_date)

    def test_get_date_of_most_recent_ena_sample(self):
        most_recent_sample_date = None
        # adds 50 loaded and 50 failed samples, the date should only take into account loaded samples
        samples = [get_mocked_ena_sample(faker=self.faker) for _ in range(50)] + \
                  [get_mocked_ena_sample(faker=self.faker, job_status=JobStatus.FAILED_LOAD) for _ in range(50)]
        for sample_ena, sample, job in samples:
            self.session.add_all([sample_ena, sample, job])
            if job.status == JobStatus.LOADED:
                if most_recent_sample_date is None:
                    most_recent_sample_date = sample_ena.first_created
                if sample_ena.first_created > most_recent_sample_date:
                    most_recent_sample_date = sample_ena.first_created
        self.session.commit()
        observed_date = get_date_of_most_recent_ena_sample(session=self.session)
        self.assertEqual(observed_date, most_recent_sample_date.date())

    def test_get_date_of_most_recent_ena_sample_empty(self):
        observed_date = get_date_of_most_recent_ena_sample(session=self.session)
        self.assertIsNone(observed_date)

    def test_get_date_of_last_ena_check(self):
        logs = [get_mocked_log(faker=self.faker) for _ in range(50)]
        self.session.add_all(logs)
        self.session.commit()
        observed_date = get_date_of_last_check(session=self.session, data_source=DataSource.ENA)
        self.assertIsNotNone(observed_date)

    def test_get_date_of_last_ena_check_empty(self):
        observed_date = get_date_of_last_check(session=self.session, data_source=DataSource.ENA)
        self.assertIsNone(observed_date)
        observed_date = get_date_of_last_check(session=self.session, data_source=DataSource.GISAID)
        self.assertIsNone(observed_date)

    def test_get_date_of_last_ena_update(self):
        # NOTE: the implementation that works in Postgres does not work in SQLite!!
        logs = [get_mocked_log(faker=self.faker, source=DataSource.ENA) for _ in range(50)]
        self.session.add_all(logs)
        self.session.commit()
        observed_date = get_date_of_last_update(session=self.session, data_source=DataSource.ENA)
        self.assertIsNotNone(observed_date)

    def test_get_date_of_last_ena_update_empty(self):
        observed_date = get_date_of_last_update(session=self.session, data_source=DataSource.ENA)
        self.assertIsNone(observed_date)
        observed_date = get_date_of_last_update(session=self.session, data_source=DataSource.GISAID)
        self.assertIsNone(observed_date)
