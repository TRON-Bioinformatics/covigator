from unittest import TestCase

from faker import Faker

from covigator.database.database import Database
from covigator.database.precomputed import Precomputer
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.faked_objects import FakeConfiguration
from covigator.tests.unit_tests.mocked import mock_samples_and_variants


class TestPrecomputer(AbstractTest):

    def setUp(self) -> None:
        mock_samples_and_variants(session=self.session, faker=self.faker, num_samples=100)
        self.precomputer = Precomputer(session=self.session)

    def test_load_dn_ds(self):
        self.precomputer.load_dn_ds()
