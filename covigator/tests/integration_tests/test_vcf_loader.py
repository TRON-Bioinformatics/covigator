from unittest import TestCase
import pkg_resources
import covigator.tests
from covigator.database.database import Database
from covigator.database.model import Variant, VariantObservation, DataSource, Sample
from covigator.pipeline.vcf_loader import VcfLoader


class VcfLoaderTests(TestCase):

    def test_vcf_loader_all_variants(self):
        session = Database(test=True).get_database_session()
        vcf_file = pkg_resources.resource_filename(covigator.tests.__name__, "resources/Sars_cov_2.all_variants.snpeff.vcf.gz")
        sample = Sample(id="TEST1", source=DataSource.ENA)
        VcfLoader().load(vcf_file, sample, session)
        session.commit()
        self.assertEqual(session.query(Variant).count(), 89709)
        self.assertEqual(session.query(VariantObservation).count(), 89709)
