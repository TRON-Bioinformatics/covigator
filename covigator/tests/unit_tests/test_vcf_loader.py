import pkg_resources
import covigator.tests
from covigator.database.model import Variant, VariantObservation, DataSource, SubclonalVariantObservation, \
    SampleEna, LowFrequencyVariantObservation, SubclonalVariant, \
    LowFrequencyVariant, LowQualityClonalVariant, LowQualityClonalVariantObservation, VariantCovid19Portal, \
    VariantObservationCovid19Portal, SampleCovid19Portal
from covigator.pipeline.vcf_loader import VcfLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest


class VcfLoaderTests(AbstractTest):

    def setUp(self) -> None:
        self.sample_ena = SampleEna(run_accession="TEST1", fastq_ftp="something", fastq_md5="else", num_fastqs=2,
                                    covered_bases=30000, read_count=100000)
        # these samples are not eligible to intrahost mutations
        self.sample_ena_low_covered_bases = SampleEna(
            run_accession="TEST2", fastq_ftp="something", fastq_md5="else", num_fastqs=2,
            covered_bases=100, read_count=100000)
        self.sample_ena_low_read_counts = SampleEna(
            run_accession="TEST3", fastq_ftp="something", fastq_md5="else", num_fastqs=2,
            covered_bases=30000, read_count=100)
        self.sample_c19dp = SampleCovid19Portal(run_accession="TEST4", fasta_url="something")
        self.session.add_all([
            self.sample_ena, self.sample_c19dp,
            self.sample_ena_low_covered_bases, self.sample_ena_low_read_counts])
        self.session.commit()

    def test_vcf_loader_ena(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()

        # check that one variant of each type was loaded
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 2)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 2)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 4)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 4)
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
        self.assertEqual(variant_observation.dp, 110)
        self.assertAlmostEqual(variant_observation.vaf, 0.2)

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

    def test_sample_low_read_counts(self):
        """
        No intrahost variants should be loaded if the sample has low read counts
        """
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena_low_covered_bases.run_accession, source=DataSource.ENA, session=self.session)
        self.session.commit()

        # check that one variant of each type was loaded
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 0)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 0)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 6)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 6)
        self.assertEqual(self.session.query(LowQualityClonalVariant).count(), 1)
        self.assertEqual(self.session.query(LowQualityClonalVariantObservation).count(), 1)

    def test_sample_low_covered_bases(self):
        """
        No intrahost mutations because the sample is excluded due to low covered bases
        """
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_ena_low_covered_bases.run_accession,
                         source=DataSource.ENA, session=self.session)
        self.session.commit()

        # check that one variant of each type was loaded
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        self.assertEqual(self.session.query(SubclonalVariant).count(), 0)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 0)
        self.assertEqual(self.session.query(LowFrequencyVariant).count(), 6)
        self.assertEqual(self.session.query(LowFrequencyVariantObservation).count(), 6)
        self.assertEqual(self.session.query(LowQualityClonalVariant).count(), 1)
        self.assertEqual(self.session.query(LowQualityClonalVariantObservation).count(), 1)

    def test_vcf_loader_c19dp(self):
        """
        All mutations are loaded as clonal when it is a C19DP sample
        """
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, run_accession=self.sample_c19dp.run_accession, source=DataSource.COVID19_PORTAL, session=self.session)
        self.session.commit()

        # check that one variant of each type was loaded
        self.assertEqual(self.session.query(VariantCovid19Portal).count(), 8)
        self.assertEqual(self.session.query(VariantObservationCovid19Portal).count(), 8)

        # check that the variant was loaded with the correct data
        variant = self.session.query(VariantCovid19Portal).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(VariantObservationCovid19Portal).first()
        self.assertEqual(variant_observation.sample, self.sample_c19dp.run_accession)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertIsNone(variant_observation.dp)
        self.assertAlmostEqual(variant_observation.vaf, 0.9)

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
