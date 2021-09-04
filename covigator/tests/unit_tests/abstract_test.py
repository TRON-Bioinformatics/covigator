from unittest import TestCase
from faker import Faker
from logzero import logger
from sqlalchemy import MetaData

from covigator.database.database import Database
from covigator.tests.unit_tests.faked_objects import FakeConfiguration


class AbstractTest(TestCase):

    @classmethod
    def setUpClass(cls):
        """On inherited classes, run our `setUp` method"""
        # Inspired via http://stackoverflow.com/questions/1323455/python-unit-test-with-base-and-sub-class/17696807#17696807
        if cls is not AbstractTest and cls.setUp is not AbstractTest.setUp:
            orig_setUp = cls.setUp

            def setUpOverride(self, *args, **kwargs):
                AbstractTest.setUp(self)
                return orig_setUp(self, *args, **kwargs)

            cls.setUp = setUpOverride

    def setUp(self) -> None:
        self.config = FakeConfiguration()
        self.database = Database(test=True, config=self.config)
        self.session = self.database.get_database_session()
        self.faker = Faker()

    def tearDown(self) -> None:
        logger.info("Cleaning the database")
        try:
            meta = MetaData()
            meta.drop_all(bind=self.database.engine)
        except Exception as e:
            logger.error("Error cleaning the database")
            logger.exception(e)
        #for table in reversed(meta.sorted_tables):
        #    self.session.execute(table.delete())
        #self.session.commit()

