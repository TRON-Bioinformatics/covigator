from unittest import TestCase

from covigator.misc.country_parser import CountryParser


class TestCountryParser(TestCase):

    def test_country_parser(self):
        parser = CountryParser()
        italy = parser.parse_country("iTAly")
        self.assertEqual(italy.country_raw, "iTAly")
        self.assertEqual(italy.country, "Italy")
        self.assertEqual(italy.country_alpha_2, "IT")
        self.assertEqual(italy.country_alpha_3, "ITA")
        self.assertEqual(italy.continent, "Europe")
        self.assertEqual(italy.continent_alpha_2, "EU")

        sofa = parser.parse_country("my sofa")
        self.assertEqual(sofa.country_raw, "my sofa")
        self.assertEqual(sofa.country, "Not available")
        self.assertEqual(sofa.country_alpha_2, "None")
        self.assertEqual(sofa.country_alpha_3, "None")
        self.assertEqual(sofa.continent, "None")
        self.assertEqual(sofa.continent_alpha_2, "None")

        empty = parser.parse_country("  ")
        self.assertEqual(empty.country_raw, "my sofa")
        self.assertEqual(empty.country, "Not available")
        self.assertEqual(empty.country_alpha_2, "None")
        self.assertEqual(empty.country_alpha_3, "None")
        self.assertEqual(empty.continent, "None")
        self.assertEqual(empty.continent_alpha_2, "None")
