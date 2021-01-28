from unittest import TestCase
import pkg_resources
import covigator.tests
from covigator.database import Database
from covigator.model import Variant, VariantObservation
from covigator.processor.vcf_loader import VcfLoader


class VcfLoaderTests(TestCase):

    def test_vcf_loader(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff.vcf")
        sample = "TEST1"
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 1)
        self.assertEqual(session.query(VariantObservation).count(), 1)
        variant = session.query(Variant).first()
        self.assertEqual(variant.chromosome, "MN908947.3")
        self.assertEqual(variant.position, 23403)
        self.assertEqual(variant.reference, "A")
        self.assertEqual(variant.alternate, "G")
        self.assertEqual(variant.gene_name, "S")
        variant_observation = session.query(VariantObservation).first()
        self.assertEqual(variant_observation.sample, sample)
        self.assertEqual(variant_observation.chromosome, "MN908947.3")
        self.assertEqual(variant_observation.position, 23403)
        self.assertEqual(variant_observation.reference, "A")
        self.assertEqual(variant_observation.alternate, "G")

    def test_vcf_loader_without_dp4(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_dp4.vcf")
        sample = "TEST1"
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 1)
        self.assertEqual(session.query(VariantObservation).count(), 1)

    def test_vcf_loader_without_annotation(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_without_annotation.vcf")
        sample = "TEST1"
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 1)
        self.assertEqual(session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_insertion(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_insertion.vcf")
        sample = "TEST1"
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 1)
        self.assertEqual(session.query(VariantObservation).count(), 1)

    def test_vcf_loader_with_deletion(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/snpeff_with_deletion.vcf")
        sample = "TEST1"
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 1)
        self.assertEqual(session.query(VariantObservation).count(), 1)
