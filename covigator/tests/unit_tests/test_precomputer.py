from sqlalchemy import and_, func

from covigator.database.model import PrecomputedSynonymousNonSynonymousCounts, RegionType, DataSource, \
    PrecomputedOccurrence
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.precomputations.load_top_occurrences import TopOccurrencesLoader
from covigator.precomputations.loader import PrecomputationsLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, MOCKED_GENES, MOCKED_DOMAINS


class TestPrecomputer(AbstractTest):

    def setUp(self) -> None:
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        self.ns_counts_loader = NsSCountsLoader(session=self.session)
        self.top_occurrences_loader = TopOccurrencesLoader(session=self.session)
        self.precomputations_loader = PrecomputationsLoader(session=self.session)

    def test_load_dn_ds(self):
        self.ns_counts_loader.load()
        for g in MOCKED_GENES:
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == g)).count(),
                0)
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == g,
                         PrecomputedSynonymousNonSynonymousCounts.source == DataSource.ENA.name)).count(),
                0)
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == g,
                         PrecomputedSynonymousNonSynonymousCounts.source == DataSource.GISAID.name)).count(),
                0)
            self.assertEqual(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type != RegionType.GENE.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == g)).count(),
                0)
        self.assertGreater(
            self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.CODING_REGION.name)).count(), 0)
        self.assertGreater(
            self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.CODING_REGION.name,
                     PrecomputedSynonymousNonSynonymousCounts.source == DataSource.ENA.name)).count(), 0)
        self.assertGreater(
            self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.CODING_REGION.name,
                     PrecomputedSynonymousNonSynonymousCounts.source == DataSource.GISAID.name)).count(), 0)

        s_genes = self.session.query(func.sum(PrecomputedSynonymousNonSynonymousCounts.s)).filter(
            PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE.name).scalar()
        s_coding_region = self.session.query(func.sum(PrecomputedSynonymousNonSynonymousCounts.s)).filter(
            PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.CODING_REGION.name).scalar()
        self.assertEqual(s_genes, s_coding_region)

        ns_genes = self.session.query(func.sum(PrecomputedSynonymousNonSynonymousCounts.ns)).filter(
            PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE.name).scalar()
        ns_coding_region = self.session.query(func.sum(PrecomputedSynonymousNonSynonymousCounts.ns)).filter(
            PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.CODING_REGION.name).scalar()
        self.assertEqual(ns_genes, ns_coding_region)

        for d in MOCKED_DOMAINS:
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.DOMAIN.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == d)).count(),
                0)
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.DOMAIN.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == d,
                         PrecomputedSynonymousNonSynonymousCounts.source == DataSource.ENA.name)).count(),
                0)
            self.assertGreater(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.DOMAIN.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == d,
                         PrecomputedSynonymousNonSynonymousCounts.source == DataSource.GISAID.name)).count(),
                0)
            self.assertEqual(
                self.session.query(PrecomputedSynonymousNonSynonymousCounts).filter(
                    and_(PrecomputedSynonymousNonSynonymousCounts.region_type != RegionType.DOMAIN.name,
                         PrecomputedSynonymousNonSynonymousCounts.region_name == d)).count(),
                0)

    def test_load_precomputed_occurrences(self):
        self.assertEqual(self.session.query(PrecomputedOccurrence).count(), 0)
        self.precomputations_loader.load_table_counts()     # table counts precomputations are needed
        self.top_occurrences_loader.load()
        self.assertGreater(self.session.query(PrecomputedOccurrence).count(), 0)
        for g in MOCKED_GENES:
            occurrences = self.session.query(PrecomputedOccurrence).filter(PrecomputedOccurrence.gene_name == g).all()
            self.assertGreater(len(occurrences), 0)
            for o in occurrences:
                self.assertGreater(o.total, 0)
                self.assertGreater(o.frequency, 0.0)
                self.assertIsNotNone(o.variant_id)
                self.assertIsNotNone(o.gene_name)
                self.assertIsNotNone(o.domain)
                self.assertIsNotNone(o.annotation)
