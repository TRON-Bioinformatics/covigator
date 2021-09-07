from sqlalchemy import and_, func

from covigator.database.model import PrecomputedDnDs, RegionType
from covigator.database.precomputed import Precomputer
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, MOCKED_GENES


class TestPrecomputer(AbstractTest):

    def setUp(self) -> None:
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        self.precomputer = Precomputer(session=self.session)

    def test_load_dn_ds(self):
        self.precomputer.load_dn_ds()
        for g in MOCKED_GENES:
            self.assertGreater(
                self.session.query(PrecomputedDnDs).filter(
                    and_(PrecomputedDnDs.region_type == RegionType.GENE.name, PrecomputedDnDs.region_name == g)).count(),
                0)
            self.assertEqual(
                self.session.query(PrecomputedDnDs).filter(
                    and_(PrecomputedDnDs.region_type != RegionType.GENE.name, PrecomputedDnDs.region_name == g)).count(),
                0)
        self.assertGreater(
            self.session.query(PrecomputedDnDs).filter(
                and_(PrecomputedDnDs.region_type == RegionType.CODING_REGION.name)).count(), 0)

        s_genes = self.session.query(func.sum(PrecomputedDnDs.s)).filter(
            PrecomputedDnDs.region_type == RegionType.GENE.name).scalar()
        s_coding_region = self.session.query(func.sum(PrecomputedDnDs.s)).filter(
            PrecomputedDnDs.region_type == RegionType.CODING_REGION.name).scalar()
        self.assertEqual(s_genes, s_coding_region)

        ns_genes = self.session.query(func.sum(PrecomputedDnDs.ns)).filter(
            PrecomputedDnDs.region_type == RegionType.GENE.name).scalar()
        ns_coding_region = self.session.query(func.sum(PrecomputedDnDs.ns)).filter(
            PrecomputedDnDs.region_type == RegionType.CODING_REGION.name).scalar()
        self.assertEqual(ns_genes, ns_coding_region)
