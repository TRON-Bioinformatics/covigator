from datetime import date
from unittest import TestCase
import random

from covigator.dashboard.figures import get_accumulated_samples_by_country
from covigator.database import Database
from covigator.model import SampleEna, Sample, DataSource, JobEna, JobStatus


class FiguresTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()

    def _get_mocked_sample(self):
        identifier = str(random.uniform(1, 1000000))
        sample_ena = SampleEna(
            run_accession=identifier,
            first_created=date.fromisoformat("2021-01-01"),
            country="Germany",
            fastq_ftp="ftp.fastq",
            fastq_md5="123456789",
            num_fastqs=1
        )
        sample = Sample(
            id=identifier,
            source=DataSource.ENA,
            ena_id=identifier
        )
        job = JobEna(
            run_accession=identifier,
            status=JobStatus.LOADED
        )
        return sample_ena, sample, job

    def test_samples_by_country(self):
        # populates the ENA samples tables
        for _ in range(100):
            sample_ena, sample, job = self._get_mocked_sample()
            self.session.add_all([sample_ena, sample, job])
        self.session.commit()
        figure = get_accumulated_samples_by_country(session=self.session)
        self.assertIsNotNone(figure)

    def test_samples_by_country_no_data(self):
        figure = get_accumulated_samples_by_country(session=self.session)
        self.assertIsNone(figure)
