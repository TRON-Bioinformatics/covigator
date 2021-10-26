from covigator.database.database import Database
from covigator.database.model import Gene, get_table_versioned_name, Variant, Conservation, Domain
import pandas as pd
from covigator.configuration import Configuration
from covigator.tests.unit_tests.abstract_test import AbstractTest


class DatabaseInitialisationTests(AbstractTest):

    def test_genes_table_initialisation(self):
        database = Database(test=True, config=self.config)
        session = database.get_database_session()
        self.assertGreater(session.query(Gene).count(), 0)
        for g in session.query(Gene).all():
            self.assertIsNotNone(g.identifier)
            self.assertIsNotNone(g.name)
            self.assertGreater(g.fraction_synonymous, 0.0)
            self.assertGreater(g.fraction_non_synonymous, 0.0)
        self.assertGreater(session.query(Domain).count(), 0)
        for d in session.query(Domain).all():
            self.assertIsNotNone(d.name)
            self.assertGreater(d.fraction_synonymous, 0.0)
            self.assertGreater(d.fraction_non_synonymous, 0.0)
            self.assertIsNotNone(d.gene_identifier)
            self.assertIsNotNone(d.gene_name)

    def test_genes_table_initialisation_not_twice(self):
        database = Database(test=True, config=self.config)
        session = database.get_database_session()
        count_genes = session.query(Gene).count()
        count_domains = session.query(Domain).count()

        # creates another connection
        database2 = Database(test=True, config=self.config)
        session2 = database2.get_database_session()

        count_genes_2 = session2.query(Gene).count()
        count_domains_2 = session2.query(Domain).count()
        self.assertEqual(count_genes, count_genes_2)
        self.assertEqual(count_domains, count_domains_2)
        self.assertGreater(count_genes, 0)
        self.assertGreater(count_domains, 0)

    def test_versioned_tables(self):
        config = Configuration()
        config.db_table_version = "_v1"
        self.assertEqual("gene_v1", get_table_versioned_name('gene', config=config))
        config.db_table_version = "_v2"
        self.assertEqual("variant_v2", get_table_versioned_name('variant', config=config))

    def test_conservation_loader(self):
        database = Database(test=True, config=self.config)
        session = database.get_database_session()
        self.assertGreater(session.query(Conservation).count(), 0)
        conservation_values = pd.read_sql(session.query(Conservation).statement, session.bind)
        self.assertEqual(conservation_values.shape[1], 6)
        self.assertGreater(conservation_values.shape[0], 1000)
        self.assertEqual(conservation_values[conservation_values.conservation.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_sarbecovirus.isna()].shape[0], 0)
        self.assertEqual(conservation_values[conservation_values.conservation_vertebrates.isna()].shape[0], 0)

    def test_no_repeated_genes(self):
        database = Database(test=True, config=self.config)
        session = database.get_database_session()
        genes = session.query(Gene).all()
        gene_names = [g.name for g in genes]
        unique_gene_names = set(gene_names)
        self.assertEqual(len(gene_names), len(unique_gene_names))
