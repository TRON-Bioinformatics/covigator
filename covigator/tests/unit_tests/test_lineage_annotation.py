from covigator.database.database import Database
from covigator.database.queries import Queries
from covigator.configuration import Configuration
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.references.lineage_annotation import LineageAnnotationsLoader


class LineageAnnotationTest(AbstractTest):
    def setUp(self) -> None:
        self.queries = Queries(session=self.session)
        self.loader = LineageAnnotationsLoader(self.session)
    def test_find_gene(self):
        gene_locations = [266, 13483]
        intergenic_locations = [265, 21556]
        for loc in gene_locations:
            x = self.loader.find_gene(loc)
            self.assertEqual(x, ('ORF1ab', 266, 21555))
        for loc in intergenic_locations:
            x = self.loader.find_gene(loc)
            self.assertIsNone(x[0])

    def test_get_hgvs_from_nuc_snp(self):
        # Test for different SNP classes
        missense_snp = {"genomic_position": 21765, "reference": "T", "alternate": "C"}
        synonymous_snp = {"genomic_position": 913, "reference": "C", "alternate": "T"}
        intergenic_snp = {"genomic_position": 240, "reference": "C", "alternate": "T"}

        x = self.loader.get_hgvs_from_nuc_snp(missense_snp)
        self.assertEqual(
            x, {'hgvs_p': 'p.I68T', 'position': 68, 'reference': 'I', 'alternate': 'T',
                'annotation': 'missense_variant', 'protein': 'S'}
        )
        x = self.loader.get_hgvs_from_nuc_snp(synonymous_snp)
        self.assertEqual(
            x, {'hgvs_p': 'p.S216S', 'position': 216, 'reference': 'S', 'alternate': 'S',
                'annotation': 'synonymous_variant', 'protein': 'ORF1ab'}
        )
        x = self.loader.get_hgvs_from_nuc_snp(intergenic_snp)
        self.assertEqual(
            x, {'hgvs_p': None, 'position': 240, 'reference': 'C', 'alternate': 'T',
                'annotation': 'intergenic_variant', 'protein': None}
        )

    def test_get_hgvs_from_nuc_deletion(self):
        # Test deletions with different effects
        # 21990:TTTA>T = p.Y144del
        normal_deletion = {"genomic_position": 21990, "length": 3}
        # 21990:TTTATTACCA>T = p.Y144_H146del
        normal_deletion_long = {"genomic_position": 21990, "length": 9}
        # 22029:AGTTCAG>A = p.F157_R158del
        normal_deletion_long_2 = {"genomic_position": 22029, "length": 6}
        # 21992:TA>T = p.Y144fs
        fs_deletion = {"genomic_position": 21992, "length": 2}
        # 21992:TATTACCACAAAAACAACAAAAGTTGGATGGAAA>T = p.Y144_S155delinsC
        insdel_deletion = {"genomic_position": 21992, "length": 33}
        # 21633:TACCCCCTGCATA>T = p.L24_Y28delinsF
        insdel_deletion_2 = {"genomic_position": 21633, "length": 12}

        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(normal_deletion),
            {'hgvs_p': 'p.Y144del', 'reference': 'Y', 'alternate': 'del',
             'position': 144, "annotation": "conservative_inframe_deletion", "protein": "S"}
        )
        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(normal_deletion_long),
            {'hgvs_p': 'p.Y144_H146del', 'reference': 'YYH', 'alternate': 'del',
             'position': 144, "annotation": "disruptive_inframe_deletion", "protein": "S"}
        )
        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(normal_deletion_long_2),
            {'hgvs_p': 'p.F157_R158del', 'reference': 'FR', 'alternate': 'del',
             'position': 157, "annotation": "disruptive_inframe_deletion", 'protein': "S"}
        )
        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(fs_deletion),
            {'hgvs_p': 'p.Y144fs', 'reference': 'Y', 'alternate': 'del',
             'position': 144, 'annotation': 'frameshift_variant', 'protein': 'S'}
        )
        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(insdel_deletion),
            {'hgvs_p': 'p.Y144_S155delinsC', 'reference': 'YYHKNNKSWMES', 'alternate': 'del',
             'position': 144, "annotation": "disruptive_inframe_deletion", "protein": "S"}
        )
        self.assertEqual(
            self.loader.get_hgvs_from_nuc_deletion(insdel_deletion_2),
            {'hgvs_p': 'p.L24_Y28delinsF', 'reference': 'LPPAY', 'alternate': 'del',
             'position': 24, "annotation": "disruptive_inframe_deletion", "protein": "S"}
        )

    def test_parse_mutation_sites(self):
        # Test that pangolin mutation definitions are parsed correctly
        parsed_site = self.loader._parse_mutation_sites("nuc:C16176T")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 1)
        self.assertIsInstance(parsed_site[0], dict)

        parsed_site = self.loader._parse_mutation_sites("spike:D614G")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 1)
        self.assertIsInstance(parsed_site[0], dict)

        # Test that grouped AA level mutation are split into separate mutations
        parsed_site = self.loader._parse_mutation_sites("n:RG203KR")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 2)
        self.assertTrue(parsed_site[0]["variant_id"] == "N:R203K")
        self.assertTrue(parsed_site[1]["variant_id"] == "N:G204R")

    def test_create_constellation_pango_mapping(self):
        fake_lineage_constellation = {
            "A": {"pangolin_lineage_list": set(["A.1"]), "parent_lineage_id": None},
            "B": {"pangolin_lineage_list": set(["B.1"]), "parent_lineage_id": "A.1"}}
        mapping = self.loader._create_constellation_pango_mapping(fake_lineage_constellation)
        self.assertIsInstance(fake_lineage_constellation, dict)
        self.assertEqual(len(mapping), 2)

    def test_find_parent_sites(self):
        fake_lineage_constellation = {
            "A": {"pangolin_lineage_list": set(["A.1"]),
                  "parent_lineage_id": None,
                  "lineage_mutations": [{"variant_id": 1}, {"variant_id": 2}]},
            "B": {"pangolin_lineage_list": set(["B.1"]),
                  "parent_lineage_id": "A.1",
                  "lineage_mutations": [{"variant_id": 3}, {"variant_id": 4}]}
        }
        b_sites = self.loader._find_parent_sites("B",
            fake_lineage_constellation,
            LineageAnnotationsLoader._create_constellation_pango_mapping(fake_lineage_constellation))
        self.assertIsNotNone(b_sites)
        self.assertEqual(len(b_sites), 4)
        # Test that parent sites are returned
        self.assertTrue({"variant_id": 1} in b_sites)
        self.assertTrue({"variant_id": 2} in b_sites)


