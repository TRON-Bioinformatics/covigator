from parameterized import parameterized
from sqlalchemy import and_, func

from covigator.database.model import PrecomputedSynonymousNonSynonymousCounts, RegionType, DataSource, \
    PrecomputedOccurrence, PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, PrecomputedIndelLength, \
    PrecomputedAnnotation, PrecomputedVariantAbundanceHistogram, PrecomputedVariantsPerLineage, VariantCooccurrence
from covigator.precomputations.load_cooccurrences import CooccurrenceMatrixLoader
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.precomputations.load_top_occurrences import TopOccurrencesLoader
from covigator.precomputations.load_variants_per_lineage import VariantsPerLineageLoader
from covigator.precomputations.loader import PrecomputationsLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.mocked import mock_samples_and_variants, MOCKED_GENES, MOCKED_DOMAINS


class TestPrecomputer(AbstractTest):

    def setUp(self) -> None:
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        self.ns_counts_loader = NsSCountsLoader(session=self.session)
        self.top_occurrences_loader = TopOccurrencesLoader(session=self.session)
        self.precomputations_loader = PrecomputationsLoader(session=self.session)
        self.precomputations_lineage = VariantsPerLineageLoader(session=self.session)

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

    def test_load_table_counts(self):
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

    def test_load_counts_variants_per_sample(self):
        self.assertEqual(self.session.query(PrecomputedVariantsPerSample).count(), 0)
        self.precomputations_loader.load_counts_variants_per_sample()
        self.assertGreater(self.session.query(PrecomputedVariantsPerSample).count(), 0)
        p: PrecomputedVariantsPerSample
        for p in self.session.query(PrecomputedVariantsPerSample).all():
            self.assertGreater(p.count, 0)
            self.assertIsNotNone(p.source)
            self.assertIsNotNone(p.variant_type)
            self.assertIsNotNone(p.number_mutations)

    def test_load_count_substitutions(self):
        self.assertEqual(self.session.query(PrecomputedSubstitutionsCounts).count(), 0)
        self.precomputations_loader.load_count_substitutions()
        self.assertGreater(self.session.query(PrecomputedSubstitutionsCounts).count(), 0)
        p: PrecomputedSubstitutionsCounts
        for p in self.session.query(PrecomputedSubstitutionsCounts).all():
            self.assertGreater(p.count, 0)
            self.assertIsNotNone(p.source)
            self.assertIsNotNone(p.variant_type)
            self.assertIsNotNone(p.reference)
            self.assertIsNotNone(p.alternate)

    def test_load_indel_length(self):
        self.assertEqual(self.session.query(PrecomputedIndelLength).count(), 0)
        self.precomputations_loader.load_indel_length()
        # TODO: mock indels
        #self.assertGreater(self.session.query(PrecomputedIndelLength).count(), 0)
        #p: PrecomputedIndelLength
        #for p in self.session.query(PrecomputedIndelLength).all():
        #    self.assertGreater(p.count, 0)
        #    self.assertIsNotNone(p.source)
        #    self.assertIsNotNone(p.gene_name)
        #    self.assertIsNotNone(p.length)

    def test_load_annotation(self):
        self.assertEqual(self.session.query(PrecomputedAnnotation).count(), 0)
        self.precomputations_loader.load_annotation()
        self.assertGreater(self.session.query(PrecomputedAnnotation).count(), 0)
        p: PrecomputedAnnotation
        for p in self.session.query(PrecomputedAnnotation).all():
            self.assertGreater(p.count, 0)
            self.assertIsNotNone(p.source)
            self.assertIsNotNone(p.annotation)

    def test_load_variant_abundance_histogram(self):
        self.assertEqual(self.session.query(PrecomputedVariantAbundanceHistogram).count(), 0)
        self.precomputations_loader.load_variant_abundance_histogram()
        self.assertGreater(self.session.query(PrecomputedVariantAbundanceHistogram).count(), 0)
        p: PrecomputedVariantAbundanceHistogram
        for p in self.session.query(PrecomputedVariantAbundanceHistogram).all():
            self.assertGreaterEqual(p.count_variant_observations, 0)
            self.assertGreaterEqual(p.count_variant_observations, p.count_unique_variants)
            self.assertIsNotNone(p.source)
            self.assertIsNotNone(p.bin_size)
            self.assertIsNotNone(p.position_bin)

    def test_load_variants_per_lineage(self):
        self.assertEqual(self.session.query(PrecomputedVariantsPerLineage).count(), 0)
        self.precomputations_lineage.load()
        self.assertGreater(self.session.query(PrecomputedVariantsPerLineage).count(), 0)
        for p in self.session.query(PrecomputedVariantsPerLineage).all():
            self.assertGreater(p.count_observations, 0)
            self.assertIsNotNone(p.lineage)
            self.assertNotEqual(p.lineage, "")
            self.assertIsNotNone(p.variant_id)
            self.assertIsNotNone(p.country)

    @parameterized.expand([DataSource.ENA.name, DataSource.GISAID.name])
    def test_load_cooccurrence_matrix(self, source):

        variant_cooccurrence_klass = self.queries.get_variant_cooccurrence_klass(source)

        self.assertEqual(self.session.query(variant_cooccurrence_klass).count(), 0)
        CooccurrenceMatrixLoader(self.session, source=source).load(maximum_length=10)
        self.assertGreater(self.session.query(variant_cooccurrence_klass).count(), 0)
        found_greater_one = False
        for p in self.session.query(variant_cooccurrence_klass).all():
            self.assertGreater(p.count, 1)  # unique observations are deleted
            found_greater_one = p.count > 1 or found_greater_one
            self.assertIsNotNone(p.variant_id_one)
            self.assertIsNotNone(p.variant_id_two)
        self.assertTrue(found_greater_one)
