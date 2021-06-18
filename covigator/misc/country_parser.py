from dataclasses import dataclass

import pycountry
import pycountry_convert


@dataclass
class ParsedCountry:
    country_raw: str
    country: str
    country_alpha_2: str
    country_alpha_3: str
    continent_alpha_2: str
    continent: str


empty_country = ParsedCountry(
    country_raw="",
    country="Not available",
    country_alpha_2="None",    # don't use NA as that is the id of Namibia
    country_alpha_3="None",
    continent_alpha_2="None",
    continent="None"
)


class CountryParser:

    def __init__(self):
        self.cache = {}

    def parse_country(self, country_raw: str) -> ParsedCountry:
        try:
            if country_raw is None or country_raw.strip() == "":
                parsed_country = empty_country
            else:
                parsed_country = self.cache.get(country_raw)
                if parsed_country is None:
                    match = pycountry.countries.search_fuzzy(country_raw)[0]
                    continent = pycountry_convert.country_alpha2_to_continent_code(match.alpha_2)
                    parsed_country = ParsedCountry(
                        country_raw=country_raw,
                        country=match.name,
                        country_alpha_2=match.alpha_2,
                        country_alpha_3=match.alpha_3,
                        continent_alpha_2=continent,
                        continent=pycountry_convert.convert_continent_code_to_continent_name(continent))
                    self.cache[country_raw] = parsed_country
        except (LookupError, AttributeError):
            parsed_country = empty_country
            parsed_country.country_raw = country_raw
        return parsed_country
