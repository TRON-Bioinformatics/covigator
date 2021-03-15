from unittest import TestCase
from covigator.database import Database
from covigator.model import Gene


class FiguresTests(TestCase):

    def setUp(self) -> None:
        # intialise database
        self.database = Database(test=True)
        self.session = self.database.get_database_session()

    def test_genes_table_initialisation(self):
        self.assertGreater(self.session.query(Gene).count(), 0)

    def test_genes_table_initialisation_not_twice(self):
        count_genes = self.session.query(Gene).count()

        # creates another connection
        database2 = Database(test=True)
        session2 = self.database.get_database_session()

        count_genes_2 = session2.query(Gene).count()
        self.assertEqual(count_genes, count_genes_2)
