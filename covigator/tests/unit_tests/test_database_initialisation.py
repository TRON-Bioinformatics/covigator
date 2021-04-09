from unittest import TestCase

from covigator import ENV_COVIGATOR_TABLE_VERSION
from covigator.database.database import Database
from covigator.database.model import Gene, get_table_versioned_name, Variant
import os


class DatabaseInitialisationTests(TestCase):

    def test_genes_table_initialisation(self):
        database = Database(test=True)
        session = database.get_database_session()
        self.assertGreater(session.query(Gene).count(), 0)

    def test_genes_table_initialisation_not_twice(self):
        database = Database(test=True)
        session = database.get_database_session()
        count_genes = session.query(Gene).count()

        # creates another connection
        database2 = Database(test=True)
        session2 = database2.get_database_session()

        count_genes_2 = session2.query(Gene).count()
        self.assertEqual(count_genes, count_genes_2)

    def test_versioned_tables(self):
        os.environ[ENV_COVIGATOR_TABLE_VERSION] = "_v1"
        self.assertEqual("gene_v1", get_table_versioned_name(Gene.__table__.name))
        os.environ[ENV_COVIGATOR_TABLE_VERSION] = "_v2"
        self.assertEqual("variant_v2", get_table_versioned_name(Variant.__table__.name))
