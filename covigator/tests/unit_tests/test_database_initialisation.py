from covigator.database.database import Database
from covigator.database.model import Gene, get_table_versioned_name, Conservation, Domain, Lineages
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
        count_lineages = session.query(Lineages).count()

        # creates another connection
        database2 = Database(test=True, config=self.config)
        session2 = database2.get_database_session()

        count_genes_2 = session2.query(Gene).count()
        count_domains_2 = session2.query(Domain).count()
        count_lineages_2 = session2.query(Lineages).count()

        self.assertEqual(count_genes, count_genes_2)
        self.assertEqual(count_domains, count_domains_2)
        self.assertEqual(count_genes, count_lineages_2)
        self.assertGreater(count_genes, 0)
        self.assertGreater(count_domains, 0)
        self.assertGreater(count_lineages, 0)

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

    def test_lineage_initialisation(self):
        database = Database(test=True, config=self.config)
        session = database.get_database_session()
        self.assertGreater(session.query(Lineages).count(), 0)
        lineages = session.query(Lineages).all()
        for d in lineages:
            self.assertIsNotNone(d.pango_lineage_id)
            self.assertIsNotNone(d.constellation_id)
            self.assertIsNotNone(d.variant_of_concern)
            self.assertIsNotNone(d.variant_under_investigation)

        # Test no repeated lineages --> constellation files for some lineages with an additional mutation
        pango_ids = [x.pango_lineage_id for x in lineages]
        unique_pango_ids = set(pango_ids)
        self.assertEqual(len(pango_ids), len(unique_pango_ids))

        parent_query = session.query(Lineages).filter(Lineages.parent_lineage_id.is_not(None))
        lineages_with_parents = pd.read_sql(parent_query.statement, session.bind)
        self.assertEqual(lineages_with_parents.shape[0], 10)
        self.assertEqual(lineages_with_parents.shape[1], 10)
        # Check that parents are all present in lineage table
        self.assertTrue(all([x in pango_ids for x in lineages_with_parents.parent_lineage_id]))

        # Check that parents were correctly parsed
        omicron = lineages_with_parents.loc[lineages_with_parents["who_label"] == "Omicron"]
        self.assertTrue(all([x in ["B.1.1.529", "XE-parent2"] for x in omicron.parent_lineage_id]))


