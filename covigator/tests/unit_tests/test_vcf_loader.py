import pkg_resources
import covigator.tests
from covigator.database.model import Variant, VariantObservation, Sample, DataSource, SubclonalVariantObservation, \
    SampleGisaid, SampleEna
from covigator.exceptions import CovigatorExcludedSampleTooManyMutations
from covigator.pipeline.vcf_loader import VcfLoader
from covigator.tests.unit_tests.abstract_test import AbstractTest


class VcfLoaderTests(AbstractTest):

    def setUp(self) -> None:
        self.sample = Sample(id="TEST1", source=DataSource.ENA)
        sample_ena = SampleEna(run_accession="TEST1", fastq_ftp="something", fastq_md5="else", num_fastqs=2)
        self.sample2 = Sample(id="TEST2", source=DataSource.GISAID)
        sample_gisaid = SampleGisaid(run_accession="TEST2")
        self.session.add_all([self.sample, self.sample2, sample_ena, sample_gisaid])
        self.session.commit()

    def test_vcf_loader(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, self.sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 3)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        variant = self.session.query(Variant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(VariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample.id)
        self.assertEqual(variant_observation.source, self.sample.source)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertEqual(variant_observation.dp, 38)
        self.assertEqual(variant_observation.vaf, 1.0)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 2)

    def test_vcf_loader_without_dp4(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_dp4.vcf")
        VcfLoader().load(vcf_file, self.sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_without_annotation(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_annotation.vcf")
        VcfLoader().load(vcf_file, self.sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_insertion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_insertion.vcf")
        VcfLoader().load(vcf_file, self.sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_deletion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_deletion.vcf")
        VcfLoader().load(vcf_file, self.sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_gisaid(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        VcfLoader().load(vcf_file, self.sample2, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 3)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)
        variant = self.session.query(Variant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = self.session.query(VariantObservation).first()
        self.assertEqual(variant_observation.sample, self.sample2.id)
        self.assertEqual(variant_observation.source, self.sample2.source)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 2)

    def test_too_many_mutations(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        self.session.commit()
        try:
            VcfLoader().load(vcf_file, self.sample2, self.session, max_snvs=0)
            self.assertTrue(False)
        except CovigatorExcludedSampleTooManyMutations:
            self.assertTrue(True)
