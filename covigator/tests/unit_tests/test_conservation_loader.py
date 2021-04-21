from unittest import TestCase
import pandas as pd
from covigator.database.database import Database
from covigator.database.model import Conservation
from covigator.references.conservation import ConservationLoader


class TestConservationLoader(TestCase):

    def setUp(self) -> None:
        self.session = Database(test=True).get_database_session()

    def test_conservation_loader(self):
        ConservationLoader(session=self.session).load_data()
        self.assertGreater(self.session.query(Conservation).count(), 0)
        conservation_values = pd.read_sql(self.session.query(Conservation).statement, self.session.bind)
        self.assertEqual(conservation_values.shape[1], 6)
        self.assertGreater(conservation_values.shape[0], 1000)
        self.assertEqual(conservation_values[conservation_values.conservation.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_sarbecovirus.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_vertebrates.isna()].shape[0], 0)
