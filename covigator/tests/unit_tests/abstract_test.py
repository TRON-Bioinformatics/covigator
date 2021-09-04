from unittest import TestCase
from faker import Faker
from sqlalchemy import MetaData

from covigator.database.database import Database
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class AbstractTest(TestCase):

    def setUp(self) -> None:
        self.config = FakeConfiguration()
        self.database = Database(test=True, config=self.config)
        self.session = self.database.get_database_session()
        self.faker = Faker()

    def tearDown(self) -> None:
        meta = MetaData()
        for table in reversed(meta.sorted_tables):
            self.session.execute(table.delete())
        self.session.commit()

