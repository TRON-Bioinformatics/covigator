from datetime import date
from unittest import TestCase

from sqlalchemy import and_

from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.database.model import SampleGisaid, Sample, JobGisaid, Log, DataSource, CovigatorModule
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID


class FakeGisaidAccessor(GisaidAccessor):

    def __init__(self, results, database=None):
        # uses an in memory database or the one provided
        super().__init__(input_file=None, host_tax_id=HOMO_SAPIENS_TAXID,
                         database=database if database else Database(test=True, config=Configuration()))
        self.results = results

    def _get_gisaid_runs(self):
        return self.results


class GisaidAccessorTests(TestCase):

    def test_filtering_by_instrument_platform(self):
        gisaid_accessor = FakeGisaidAccessor([
            {"run_accession": "EPI_ISL_417090",
             "virus_name": "hCoV-19/USA/WA-S37/2020",
             "instrument_platform": "illumina",
             "host_tax_id": "9606"},
            {"run_accession": "EPI_ISL_413513",
             "virus_name": "hCoV-19/South Korea/KUMC03/2020",
             "instrument_platform": "illumina",
             "host_tax_id": "9606"},
            {"run_accession": "EPI_ISL_413594",
             "virus_name": "hCoV-19/Australia/NSW08/2020",
             "instrument_platform": "ion torrent",
             "host_tax_id": "9606"},
             {"run_accession": "EPI_ISL_413594",
             "virus_name": "hCoV-19/Australia/NSW08/2020",
             "instrument_platform": "none",
             "host_tax_id": "9606"}
        ])
        gisaid_accessor.access()
        self.assertEqual(gisaid_accessor.included, 3)
        self.assertEqual(gisaid_accessor.excluded, 1)
        #self.assertEqual(gisaid_accessor.excluded_samples_by_instrument_platform.get("illumina"), 1)
        self.assertEqual(gisaid_accessor.excluded_samples_by_instrument_platform.get("none"), 1)

    def test_get_gisaid_runs(self):
        gisaid_accessor = GisaidAccessor(input_file=None, host_tax_id=HOMO_SAPIENS_TAXID,
                         database=Database(test=True, config=Configuration()))

        gisaid_accessor.access()