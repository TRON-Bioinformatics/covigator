from datetime import date, datetime
from typing import Union
import requests
from covigator.misc.country_parser import CountryParser
from covigator.database.model import SampleEna
from covigator.accessor import MINIMUM_DATE
from covigator.exceptions import CovigatorExcludedSampleTooEarlyDateException
from covigator.misc import backoff_retrier

NUMBER_RETRIES = -1


class SampleCovid19:
    pass


class AbstractAccessor:

    def __init__(self):
        self.country_parser = CountryParser()

        # this ensures there is a retry mechanism in place with a limited number of retries
        self.session = requests.session()
        self.get_with_retries = backoff_retrier.wrapper(self.session.get, NUMBER_RETRIES)

    def _parse_country(self, sample: Union[SampleEna, SampleCovid19]):
        parsed_country = self.country_parser.parse_country(
            sample.country.split(":")[0] if sample.country else "")
        sample.country_raw = sample.country
        sample.country = parsed_country.country
        sample.country_alpha_2 = parsed_country.country_alpha_2
        sample.country_alpha_3 = parsed_country.country_alpha_3
        sample.continent_alpha_2 = parsed_country.continent_alpha_2
        sample.continent = parsed_country.continent

    def _parse_dates(self, sample: Union[SampleEna, SampleCovid19]):
        sample.collection_date = self._parse_abstract(sample.collection_date, date.fromisoformat)
        sample.first_created = self._parse_abstract(sample.first_created, date.fromisoformat)
        if sample.collection_date is not None and sample.collection_date < MINIMUM_DATE:
            raise CovigatorExcludedSampleTooEarlyDateException

    def _parse_abstract(self, value, type):
        try:
            value = type(value)
        except (ValueError, TypeError):
            value = None
        return value
