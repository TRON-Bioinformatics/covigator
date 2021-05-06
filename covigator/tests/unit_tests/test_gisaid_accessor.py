from datetime import date
from unittest import TestCase
import pkg_resources
import covigator
from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.database.database import Database
from covigator.database.model import SampleGisaid, Sample, JobGisaid
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class GisaidAccessorTests(TestCase):

    def setUp(self) -> None:
        self.config = FakeConfiguration()
        self.input_fasta = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/gisaid_allprot1208.fasta")
        self.database = Database(test=True, config=self.config)

    def test_gisaid_accessor(self):
        GisaidAccessor(input_file=self.input_fasta, database=self.database).access()
        session = self.database.get_database_session()
        count_samples_gisaid = session.query(SampleGisaid).count()
        self.assertEqual(count_samples_gisaid, 2)
        count_samples = session.query(Sample).count()
        self.assertEqual(count_samples, 2)
        count_jobs = session.query(JobGisaid).count()
        self.assertEqual(count_jobs, 2)
