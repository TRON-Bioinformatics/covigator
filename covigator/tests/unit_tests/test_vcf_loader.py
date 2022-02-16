import pkg_resources
import covigator.tests
from covigator.database.model import Variant, VariantObservation, DataSource, SubclonalVariantObservation, \
    SampleGisaid, SampleEna, GisaidVariant, GisaidVariantObservation, LowFrequencyVariantObservation, SubclonalVariant, \
    LowFrequencyVariant
from covigator.pipeline.vcf_loader import VcfLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest


class VcfLoaderTests(AbstractTest):

    def setUp(self) -> None:
        self.sample_ena = SampleEna(run_accession="TEST1", fastq_ftp="something", fastq_md5="else", num_fastqs=2)
        self.sample_gisaid = SampleGisaid(run_accession="TEST2")
        self.session.add_all([self.sample_ena, self.sample_gisaid])
        self.session.commit()

    def test_vcf_loader(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 1)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 1)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 1)
        variant = self.session.query(Variant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(VariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample_ena.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertIsNone(variant_observation.dp)
        self.assertIsNone(variant_observation.vaf)

    def test_vcf_loader_without_dp4(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_dp4.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_without_annotation(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_annotation.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_insertion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_insertion.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_deletion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_deletion.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_gisaid(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_gisaid.run_accession, source=DataSource.GISAID, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(GisaidVariant).count(), 3)
        self.assertEqual(self.session.query(GisaidVariantObservation).count(), 3)
        variant = self.session.query(GisaidVariant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(GisaidVariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample_gisaid.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")

    def test_vcf_loader_vafator(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.lofreq.vcf.gz")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 2)
        self.assertEqual(self.session.query(VariantObservation).count(), 2)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 1)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 8)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 8)

        # NOTE: this test VCF was computed with thresholds of 0.2 and 0.8
        for vo in self.session.query(VariantObservation).all():
            self.assertGreaterEqual(vo.vaf, 0.8)
            self.assertGreater(vo.dp, 0)
            self.assertGreater(vo.ac, 0)

        for vo in self.session.query(SubclonalVariantObservation).all():
            self.assertGreaterEqual(vo.vaf, 0.2)
            self.assertLess(vo.vaf, 0.8)
            self.assertGreater(vo.dp, 0)
            self.assertGreater(vo.ac, 0)

        for vo in self.session.query(LowFrequencyVariantObservation).all():
            self.assertLess(vo.vaf, 0.2)
            self.assertGreater(vo.dp, 0)
            self.assertGreater(vo.ac, 0)

    def test_vcf_loader_gisaid_2(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/test.assembly.vcf.gz")
        VcfLoader().load(vcf_file, run_accession=self.sample_gisaid.run_accession, source=DataSource.GISAID, session=self.session)
        self.session.commit()
        self.assertEqual(self.session.query(GisaidVariant).count(), 13)
        self.assertEqual(self.session.query(GisaidVariantObservation).count(), 13)

        # NOTE: this test VCF was computed with thresholds of 0.2 and 0.8
        for vo in self.session.query(GisaidVariantObservation).all():
            self.assertIsNone(vo.vaf)
            self.assertIsNone(vo.dp)
            self.assertIsNone(vo.ac)
