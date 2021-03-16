from unittest import TestCase
from covigator.dashboard.figures import get_accumulated_samples_by_country, get_variants_plot
from covigator.database.database import Database
from faker import Faker

from covigator.tests.unit_tests.mocked import get_mocked_ena_sample


class FiguresTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()
        self.faker = Faker()

    def test_samples_by_country(self):
        # populates the ENA samples tables
        for _ in range(100):
            sample_ena, sample, job = get_mocked_ena_sample(faker=self.faker)
            self.session.add_all([sample_ena, sample, job])
        self.session.commit()
        figure = get_accumulated_samples_by_country(session=self.session)
        self.assertIsNotNone(figure)

    def test_samples_by_country_no_data(self):
        figure = get_accumulated_samples_by_country(session=self.session)
        self.assertIsNone(figure)

    def test_needle_plot_no_data(self):
        figure = get_variants_plot(session=self.session)
        self.assertIsNone(figure)
