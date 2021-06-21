import os
from unittest import TestCase
import pkg_resources

import covigator
from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.database.database import Database
from covigator.database.model import SampleGisaid, Sample, JobGisaid
from covigator.database.queries import Queries
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class GisaidProcessorTests(TestCase):

    def setUp(self) -> None:
        self.config = FakeConfiguration()
        self.input_fasta = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/gisaid_allprot1208.fasta")
        self.input_fasta_2 = pkg_resources.resource_filename(
            covigator.tests.__name__, "resources/gisaid_allprot1208.with_mutations.fasta")
        self.database = Database(test=True, config=self.config, initialize=True)

    def test_gisaid_pipeline_no_calls(self):
        GisaidAccessor(input_file=self.input_fasta, database=self.database).access()
        session = self.database.get_database_session()
        self.assertEqual(session.query(SampleGisaid).count(), 2)
        self.assertEqual(session.query(Sample).count(), 2)
        self.assertEqual(session.query(JobGisaid).count(), 2)

        vcf_file = GisaidPipeline(config=self.config).run(sample=session.query(SampleGisaid).first())
        self.assertTrue(os.path.exists(vcf_file))

    def test_gisaid_pipeline_some_calls(self):
        GisaidAccessor(input_file=self.input_fasta_2, database=self.database).access()
        session = self.database.get_database_session()
        self.assertEqual(session.query(SampleGisaid).count(), 2)
        self.assertEqual(session.query(Sample).count(), 2)
        self.assertEqual(session.query(JobGisaid).count(), 2)

        vcf_file = GisaidPipeline(config=self.config).run(sample=session.query(SampleGisaid).first())
        self.assertTrue(os.path.exists(vcf_file))
