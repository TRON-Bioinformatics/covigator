from unittest import TestCase

from covigator.configuration import Configuration
from covigator.references.gene_annotations import GeneAnnotationsLoader


class GeneAnnotationsLoaderTest(TestCase):

    def setUp(self) -> None:
        self.loader = GeneAnnotationsLoader()

    def test_load(self):
        self.loader.load_data()
