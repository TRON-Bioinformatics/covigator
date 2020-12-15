from unittest import TestCase
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.model import Database


class EnaDownloaderTests(TestCase):

    def test_access(self):
        sarscov2_tax_id = "2697049"
        homo_sapiens_tax_id = "9606"
        accessor = EnaAccessor(tax_id=sarscov2_tax_id, host_tax_id=homo_sapiens_tax_id, database=Database())
        accessor.access()
