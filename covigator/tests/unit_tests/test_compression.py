from unittest import TestCase
from covigator.misc.compression import compress_sequence, decompress_sequence


class TestCompression(TestCase):

    def test_compression(self):
        test_sequence = "ACGTACGTACGT"
        compressed_sequence = compress_sequence(test_sequence)
        decompressed_sequence = decompress_sequence(compressed_sequence)
        self.assertEqual(test_sequence, decompressed_sequence)