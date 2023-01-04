import pkg_resources
import covigator.tests
from covigator.database.model import Variant, VariantObservation, DataSource, SubclonalVariantObservation, \
    SampleEna, LowFrequencyVariantObservation, SubclonalVariant, \
    LowFrequencyVariant, LowQualityClonalVariant, LowQualityClonalVariantObservation
from covigator.pipeline.vcf_loader import VcfLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest


class VcfLoaderTests(AbstractTest):

    def setUp(self) -> None:
        self.sample_ena = SampleEna(run_accession="TEST1", fastq_ftp="something", fastq_md5="else", num_fastqs=2)
        self.session.add_all([self.sample_ena])
        self.session.commit()

    def test_vcf_loader(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()

        # check that one variant of each type was loaded
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 1)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 1)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 1)
        self.assertEqual(self.session.query(LowQualityClonalVariant).count(), 1)
        self.assertEqual(self.session.query(LowQualityClonalVariantObservation).count(), 1)

        # check that the variant was loaded with the correct data
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
        self.assertAlmostEqual(variant_observation.vaf, 0.9)

        variant = self.session.query(SubclonalVariant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "C")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(SubclonalVariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample_ena.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "C")
        self.assertIsNone(variant_observation.dp)
        self.assertAlmostEqual(variant_observation.vaf, 0.4)

        variant = self.session.query(LowFrequencyVariant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "T")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(LowFrequencyVariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample_ena.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "T")
        self.assertIsNone(variant_observation.dp)
        self.assertAlmostEqual(variant_observation.vaf, 0.01)

        variant = self.session.query(LowQualityClonalVariant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "ACCC")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(LowQualityClonalVariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample_ena.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "ACCC")
        self.assertIsNone(variant_observation.dp)
        self.assertAlmostEqual(variant_observation.vaf, 0.6)

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
