from sqlalchemy import and_, func

from covigator.database.model import PrecomputedSynonymousNonSynonymousCounts, RegionType, DataSource
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, MOCKED_GENES, MOCKED_DOMAINS


class TestPrecomputer(AbstractTest):

    def setUp(self) -> None:
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        self.loader = NsSCountsLoader(session=self.session)

    def test_load_dn_ds(self):
        self.loader.load()
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
