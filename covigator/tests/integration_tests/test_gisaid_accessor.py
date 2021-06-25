from datetime import date
import os
from unittest import TestCase
import pkg_resources
from Bio import SeqIO

import covigator
from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.database.model import SampleGisaid, Sample, JobGisaid
from covigator.misc.compression import decompress_sequence
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class GisaidAccessorTests(TestCase):

    def setUp(self) -> None:
        self.config = FakeConfiguration()
        dirname=os.path.dirname
        self.input_fasta = os.path.join(dirname(dirname(__file__)), "resources", "gisaid_allprot0526.fasta")
        self.input_metadata = os.path.join(dirname(dirname(__file__)), "resources", "metadata_2021_05_26.csv")
        self.database = Database(test=True, config=self.config)

    def test_gisaid_accessor(self):
        GisaidAccessor(input_fasta=self.input_fasta, input_metadata=self.input_metadata, database=self.database).access()
        session = self.database.get_database_session()
        self.assertEqual(session.query(SampleGisaid).count(), 1)
        self.assertEqual(session.query(Sample).count(), 1)
        self.assertEqual(session.query(JobGisaid).count(), 1)
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
                decompressed_sequence = decompress_sequence(seq)
                found = False
                for record in SeqIO.parse(self.input_fasta, "fasta"):
                    if str(record.seq) == decompressed_sequence:
                        found = True
                self.assertTrue(found)
