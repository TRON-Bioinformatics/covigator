from unittest import TestCase
from faker import Faker
from logzero import logger
from covigator.database.database import Database
from covigator.database.model import VariantCooccurrence, SampleEna, VariantObservation, \
    SampleGisaid, SubclonalVariantObservation, Variant, PrecomputedVariantAbundanceHistogram, PrecomputedTableCounts, \
    PrecomputedSynonymousNonSynonymousCounts, PrecomputedOccurrence, PrecomputedAnnotation, PrecomputedIndelLength, \
    PrecomputedSubstitutionsCounts, PrecomputedVariantsPerSample, Log, GisaidVariant, GisaidVariantObservation, \
    LowFrequencyVariantObservation, LowFrequencyVariant, SubclonalVariant, PrecomputedVariantsPerLineage, \
    GisaidVariantCooccurrence
from covigator.database.queries import Queries
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
        self.queries = Queries(session=self.session)
        self.faker = Faker()
        self._clean_test_database()

    def tearDown(self) -> None:
        logger.info("Cleaning the database")
        self._clean_test_database()

    def _clean_test_database(self):
        try:
            self._clean_table(GisaidVariantObservation)
            self._clean_table(VariantObservation)
            self._clean_table(SubclonalVariantObservation)
            self._clean_table(LowFrequencyVariantObservation)
            self._clean_table(SampleEna)
            self._clean_table(SampleGisaid)
            self._clean_table(VariantCooccurrence)
            self._clean_table(GisaidVariantCooccurrence)
            self._clean_table(Variant)
            self._clean_table(SubclonalVariant)
            self._clean_table(LowFrequencyVariant)
            self._clean_table(GisaidVariant)
            self._clean_table(PrecomputedVariantsPerSample)
            self._clean_table(PrecomputedSubstitutionsCounts)
            self._clean_table(PrecomputedIndelLength)
            self._clean_table(PrecomputedAnnotation)
            self._clean_table(PrecomputedOccurrence)
            self._clean_table(PrecomputedSynonymousNonSynonymousCounts)
            self._clean_table(PrecomputedTableCounts)
            self._clean_table(PrecomputedVariantAbundanceHistogram)
            self._clean_table(PrecomputedVariantsPerLineage)
            self._clean_table(Log)
        except Exception as e:
            logger.error("Error cleaning the database")
            logger.exception(e)

    def _clean_table(self, clazz):
        self.session.query(clazz).delete()
        self.session.commit()

