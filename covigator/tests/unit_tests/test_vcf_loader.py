from unittest import TestCase
import pkg_resources
import covigator.tests
from covigator.database.database import Database
from covigator.database.model import Variant, VariantObservation, Sample, DataSource, SubclonalVariantObservation
from covigator.pipeline.vcf_loader import VcfLoader
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class VcfLoaderTests(TestCase):

    def setUp(self) -> None:
        self.session = Database(test=True, config=FakeConfiguration()).get_database_session()

    def test_vcf_loader(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, self.session)
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
        self.assertEqual(variant_observation.sample, sample.id)
        self.assertEqual(variant_observation.source, sample.source)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertEqual(variant_observation.dp, 38)
        self.assertEqual(variant_observation.vaf, 1.0)
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 2)

    def test_vcf_loader_without_dp4(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_dp4.vcf")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_without_annotation(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_annotation.vcf")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_insertion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_insertion.vcf")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_deletion(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_deletion.vcf")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, self.session)
        self.session.commit()
        self.assertEqual(self.session.query(Variant).count(), 1)
        self.assertEqual(self.session.query(VariantObservation).count(), 1)

    def test_vcf_loader_gisaid(self):
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        sample = Sample(id="TEST1", source=DataSource.GISAID)
        VcfLoader().load(vcf_file, sample, self.session)
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
        self.assertEqual(variant_observation.sample, sample.id)
        self.assertEqual(variant_observation.source, sample.source)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")
        self.assertEqual(self.session.query(SubclonalVariantObservation).count(), 2)
