from unittest import TestCase
from covigator.database.database import Database
from covigator.database.model import Gene, get_table_versioned_name, Variant, Conservation
import pandas as pd
from covigator.configuration import Configuration


class DatabaseInitialisationTests(TestCase):

    def test_genes_table_initialisation(self):
        database = Database(test=True, config=Configuration())
        session = database.get_database_session()
        self.assertGreater(session.query(Gene).count(), 0)

    def test_genes_table_initialisation_not_twice(self):
        database = Database(test=True, config=Configuration())
        session = database.get_database_session()
        count_genes = session.query(Gene).count()

        # creates another connection
        database2 = Database(test=True, config=Configuration())
        session2 = database2.get_database_session()

        count_genes_2 = session2.query(Gene).count()
        self.assertEqual(count_genes, count_genes_2)

    def test_versioned_tables(self):
        config = Configuration()
        config.db_table_version = "_v1"
        self.assertEqual("gene_v1", get_table_versioned_name(Gene.__table__.name, config=config))
        config.db_table_version = "_v2"
        self.assertEqual("variant_v2", get_table_versioned_name(Variant.__table__.name, config=config))

    def test_conservation_loader(self):
        database = Database(test=True, config=Configuration())
        session = database.get_database_session()
        self.assertGreater(session.query(Conservation).count(), 0)
        conservation_values = pd.read_sql(session.query(Conservation).statement, session.bind)
        self.assertEqual(conservation_values.shape[1], 6)
        self.assertGreater(conservation_values.shape[0], 1000)
        self.assertEqual(conservation_values[conservation_values.conservation.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_sarbecovirus.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_vertebrates.isna()].shape[0], 0)
