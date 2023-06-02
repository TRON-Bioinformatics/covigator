from covigator.database.database import Database
from covigator.database.queries import Queries
from covigator.configuration import Configuration
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.references.lineage_annotation import LineageAnnotationsLoader, Hgvs, MutationInfo


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
        self.assertTrue(x.hgvs_p == 'p.I68T')
        self.assertTrue(x.position == 68)
        self.assertTrue(x.reference == "I")
        self.assertTrue(x.alternate == "T")
        self.assertTrue(x.annotation == "missense_variant")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_snp(synonymous_snp)
        self.assertTrue(x.hgvs_p == 'p.S216S')
        self.assertTrue(x.position == 216)
        self.assertTrue(x.reference == "S")
        self.assertTrue(x.alternate == "S")
        self.assertTrue(x.annotation == "synonymous_variant")
        self.assertTrue(x.protein == "ORF1ab")

        x = self.loader.get_hgvs_from_nuc_snp(intergenic_snp)
        self.assertTrue(x.hgvs_p is None)
        self.assertTrue(x.position == 240)
        self.assertTrue(x.reference == "C")
        self.assertTrue(x.alternate == "T")
        self.assertTrue(x.annotation == "intergenic_variant")
        self.assertTrue(x.protein is None)

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

        x = self.loader.get_hgvs_from_nuc_deletion(normal_deletion)
        self.assertTrue(x.hgvs_p == "p.Y144del")
        self.assertTrue(x.reference == "Y")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 144)
        self.assertTrue(x.annotation == "conservative_inframe_deletion")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_deletion(normal_deletion_long)
        self.assertTrue(x.hgvs_p == "p.Y144_H146del")
        self.assertTrue(x.reference == "YYH")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 144)
        self.assertTrue(x.annotation == "disruptive_inframe_deletion")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_deletion(normal_deletion_long_2)
        self.assertTrue(x.hgvs_p == "p.F157_R158del")
        self.assertTrue(x.reference == "FR")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 157)
        self.assertTrue(x.annotation == "disruptive_inframe_deletion")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_deletion(fs_deletion)
        self.assertTrue(x.hgvs_p == "p.Y144fs")
        self.assertTrue(x.reference == "Y")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 144)
        self.assertTrue(x.annotation == "frameshift_variant")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_deletion(insdel_deletion)
        self.assertTrue(x.hgvs_p == "p.Y144_S155delinsC")
        self.assertTrue(x.reference == "YYHKNNKSWMES")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 144)
        self.assertTrue(x.annotation == "disruptive_inframe_deletion")
        self.assertTrue(x.protein == "S")

        x = self.loader.get_hgvs_from_nuc_deletion(insdel_deletion_2)
        self.assertTrue(x.hgvs_p == "p.L24_Y28delinsF")
        self.assertTrue(x.reference == "LPPAY")
        self.assertTrue(x.alternate == "del")
        self.assertTrue(x.position == 24)
        self.assertTrue(x.annotation == "disruptive_inframe_deletion")
        self.assertTrue(x.protein == "S")


    def test_parse_mutation_sites(self):
        # Test that pangolin mutation definitions are parsed correctly
        parsed_site = self.loader._parse_mutation_sites("nuc:C16176T")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 1)
        self.assertIsInstance(parsed_site[0], MutationInfo)

        parsed_site = self.loader._parse_mutation_sites("spike:D614G")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 1)
        self.assertIsInstance(parsed_site[0], MutationInfo)

        # Test that grouped AA level mutation are split into separate mutations
        parsed_site = self.loader._parse_mutation_sites("n:RG203KR")
        self.assertIsNotNone(parsed_site)
        self.assertEqual(len(parsed_site), 2)
        self.assertTrue(parsed_site[0].variant_id == "N:R203K")
        self.assertTrue(parsed_site[1].variant_id == "N:G204R")

    def test_create_constellation_pango_mapping(self):
        fake_lineage_constellation = {
            "A": {"pangolin_lineage_list": set(["A.1"]), "parent_lineage_id": None},
            "B": {"pangolin_lineage_list": set(["B.1"]), "parent_lineage_id": "A.1"}}
        mapping = self.loader._create_constellation_pango_mapping(fake_lineage_constellation)
        self.assertIsInstance(fake_lineage_constellation, dict)
        self.assertEqual(len(mapping), 2)

    def test_find_parent_sites(self):
        fake_lineage_a_mutations = [Hgvs(hgvs_p=1, reference="A", alternate="T", position=1, annotation="A", protein="S"),
                                    Hgvs(hgvs_p=2, reference="A", alternate="T", position=1, annotation="A", protein="S")]
        fake_lineage_b_mutations = [Hgvs(hgvs_p=3, reference="A", alternate="T", position=1, annotation="A", protein="S"),
                                    Hgvs(hgvs_p=4, reference="A", alternate="T", position=1, annotation="A", protein="S")]

        fake_lineage_constellation = {
            "A": {"pangolin_lineage_list": set(["A.1"]),
                  "parent_lineage_id": None,
                  "lineage_mutations": fake_lineage_a_mutations},
            "B": {"pangolin_lineage_list": set(["B.1"]),
                  "parent_lineage_id": "A.1",
                  "lineage_mutations": fake_lineage_b_mutations}
        }
        b_sites = self.loader._find_parent_sites("B",
            fake_lineage_constellation,
            LineageAnnotationsLoader._create_constellation_pango_mapping(fake_lineage_constellation))
        self.assertIsNotNone(b_sites)
        self.assertEqual(len(b_sites), 4)
        # Test that parent sites are returned
        self.assertTrue(fake_lineage_a_mutations[0] in b_sites)
        self.assertTrue(fake_lineage_a_mutations[1] in b_sites)

    def test_is_frameshift(self):
        self.assertTrue(LineageAnnotationsLoader._is_frameshift(2))
        self.assertFalse(LineageAnnotationsLoader._is_frameshift(3))

    def test_get_protein_position(self):
        self.assertTrue(LineageAnnotationsLoader._get_protein_position(1) == 1)
        self.assertTrue(LineageAnnotationsLoader._get_protein_position(2) == 1)
        self.assertTrue(LineageAnnotationsLoader._get_protein_position(4) == 2)
        self.assertTrue(LineageAnnotationsLoader._get_protein_position(10) == 4)

    def test_get_last_position_of_deletion(self):
        self.assertTrue(LineageAnnotationsLoader._get_last_position_of_deletion(100, 2) == 100)
        self.assertTrue(LineageAnnotationsLoader._get_last_position_of_deletion(100, 3) == 101)
        self.assertTrue(LineageAnnotationsLoader._get_last_position_of_deletion(100, 4) == 100)
        self.assertTrue(LineageAnnotationsLoader._get_last_position_of_deletion(100, 6) == 102)
        self.assertTrue(LineageAnnotationsLoader._get_last_position_of_deletion(100, 9) == 103)

    def test_get_frameshift_hgvs(self):
        hgvs = LineageAnnotationsLoader._get_frameshift_hgvs("A", 1)
        self.assertIsNotNone(hgvs)
        self.assertTrue(hgvs == "p.1Afs")

    def test_get_insdel_hgvs(self):
        hgvs = LineageAnnotationsLoader._get_insdel_hgvs("A", 1, "V", 2, "G")
        self.assertIsNotNone(hgvs)
        self.assertTrue(hgvs == "p.A1_V2delinsG")

    def test_get_deletion_hgvs(self):
        hgvs = LineageAnnotationsLoader._get_deletion_hgvs("A", 1, "V", 2)
        self.assertIsNotNone(hgvs)
        self.assertTrue(hgvs == "p.A1_V2del")
        hgvs = LineageAnnotationsLoader._get_deletion_hgvs("A", 1, None, 1)
        self.assertIsNotNone(hgvs)
        self.assertTrue(hgvs == "p.A1del")

    def test_generate_variant_id(self):
        # Test nucleotide level variant id
        var_id = LineageAnnotationsLoader._generate_variant_id(None, 1, "A", "G")
        self.assertIsNotNone(var_id)
        self.assertTrue(var_id == "1:A>G")
        # Test proteomic variant id
        var_id = LineageAnnotationsLoader._generate_variant_id("S", 144, "G", "R")
        self.assertIsNotNone(var_id)
        self.assertTrue(var_id == "S:G144R")




