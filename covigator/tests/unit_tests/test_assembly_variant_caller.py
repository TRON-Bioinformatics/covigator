from unittest import TestCase

from covigator.misc.country_parser import CountryParser
from covigator.pipeline.gisaid_pipeline import AssemblyVariantCaller


class TestCountryParser(TestCase):

    def test_assembly_variant_caller(self):
        caller = AssemblyVariantCaller()
        # no mutations
        variants = caller.call_variants(sequence="ACGTACGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)
        # SNV
        variants = caller.call_variants(sequence="ACGTCCGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "A")
        self.assertEqual(snv.alternate, "C")
        self.assertEqual(snv.position, 4)
        # deletion
        variants = caller.call_variants(
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG",
            sequence="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGCAAAGTGACG")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "TGGA")
        self.assertEqual(snv.alternate, "T")
        self.assertEqual(snv.position, 45)
        # insertion
        variants = caller.call_variants(
            sequence= "CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGCCCGAGGCAAAGTGACG",
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GCCC")
        self.assertEqual(snv.position, 46)
        # another insertion
        variants = caller.call_variants(
            sequence= "CTGGTGTGAGTCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGCCCGAGGCAAAGTGACG",
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG")
        self.assertEqual(len(variants), 2)
        snv = variants[1]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GCCC")
        self.assertEqual(snv.position, 46)
        snv = variants[0]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GT")
        self.assertEqual(snv.position, 10)

