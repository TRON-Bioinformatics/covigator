from datetime import date
from unittest import TestCase
import pkg_resources
from Bio import SeqIO

import covigator
from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.database.database import Database
from covigator.database.model import SampleGisaid, Sample, JobGisaid
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
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
        self.assertEqual(session.query(SampleGisaid).count(), 2)
        self.assertEqual(session.query(Sample).count(), 2)
        self.assertEqual(session.query(JobGisaid).count(), 2)
        samples = session.query(SampleGisaid).all()
        for s in samples:
            self.assertEqual(s.country, "China")
            self.assertEqual(s.country_alpha_2, "CN")
            self.assertEqual(s.country_alpha_3, "CHN")
            self.assertEqual(s.continent, "Asia")
            self.assertEqual(s.continent_alpha_2, "AS")
            self.assertEqual(s.host, "human")
            # checks that every sequence after decompression can be found in the fasta file
            for g, seq in s.sequence.items():
                decompressed_sequence = GisaidPipeline.decompress_sequence(seq)
                found = False
                for record in SeqIO.parse(self.input_fasta, "fasta"):
                    if str(record.seq) == decompressed_sequence:
                        found = True
                self.assertTrue(found)
