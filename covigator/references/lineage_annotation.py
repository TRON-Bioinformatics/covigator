import os
import json
import re
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from sqlalchemy.orm import Session
from covigator.database.model import Lineages, LineageDefiningVariants, LineageVariant, Gene
from logzero import logger
from datetime import datetime


class LineageAnnotationsLoader:

    LINEAGE_CONSTELLATION_DIRECTORY = "constellations/constellations/definitions"
    LINEAGE_CONSTELLATION_GENOME = "constellations/constellations/data/SARS-CoV-2.json"

    def __init__(self, session: Session):
        self.session = session

        self.lineage_constellation_directory = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.LINEAGE_CONSTELLATION_DIRECTORY)
        self.gene_df = pd.read_sql(self.session.query(Gene).order_by(Gene.start).statement, self.session.bind)
        self.genome = self.load_genome()

    def load_genome(self):
        genome_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), self.LINEAGE_CONSTELLATION_GENOME)
        with open(genome_file, "r") as file_handle:
            constellation_genome = json.load(file_handle)
        return constellation_genome.get("genome")

    @staticmethod
    def _match_protein_name(protein_name: str):
        """
        This method provides a mapping between canonical protein/gene names and common synonyms.
        """
        protein_mapping = {
            "S": ["S", "s", "spike"],
            "N": ["N", "n"],
            "M": ["M", "m"],
            "E": ["E", "e"],
            "ORF1ab": ["ORF1ab", "1ab", "orf1ab", "ORF1a", "ORF1b", "orf1a", "orf1b", "nsp3", "nsp6", "nsp5", "nsp4",
                       "nsp15", "nsp13", "nsp7", "nsp12", "NSP2"],
            "ORF3a": ["ORF3a", "orf3a"],
            "ORF8": ["ORF8", "8"],
            "ORF6": ["ORF6", "orf6"]
        }
        if protein_name is None:
            return None
        for canonical_name, synonyms in protein_mapping.items():
            if protein_name in synonyms:
                return canonical_name
        return protein_name

    def find_gene(self, position):
        if position is None:
            return
        gene_coordinates = self.gene_df.query("start <= @position & end > @genomic_position").get(
            ["name", "start", "end"])
        if gene_coordinates.shape[0] == 0:
            logger.info("Position does not overlap any gene...")
            return
        gene_name = gene_coordinates.get("name").item()
        gene_name = self._match_protein_name(gene_name)
        return gene_name

    def get_hgvs_from_nuc_deletion(self, variant: dict):
        """
        Take pangolin nucleotide description and convert to hgvs notation

        ATGGT|TTA|TGAGTGATGA
        del:21990:3 --> 21990:TTTA>T -> p.Y144-

        ||||||||||||||||||||||||
          -  -  -   -   -  -  -
        GTCAATGGT|TTA|TGAGTGATG
        cds_start = 3
        deletion = del:11:3
        genomic_position = 11 - (3-1) = 9

        aa_position = 9//3

        Codon of AA position 3
        [7,8,9]
        Last base of codon --> Next protein position affected
        AA postion + 1 = 4
        Codon of position 4
        [10,11,12]
        p.hgvs of deletion = p.Y4del

        """
        genomic_position = variant.get("genomic_position")
        deletion_length = variant.get("length")
        gene_coordinates = self.gene_df.query("start <= @genomic_position & end >= @genomic_position").get(
            ["name", "start", "end"])
        if gene_coordinates.shape[0] == 0:
            logger.info("Deletion does not overlap any protein...")
            return
        protein_start = gene_coordinates.get("start").item()
        # Given a nucleotide position calculate position relative to reference protein coordinates
        # minus as positions given in df are inclusive on both sides
        position_in_cds = genomic_position - (protein_start - 1)
        # Positions are non-inclusive ->  the next genomic position is affected by the deletion
        # Get protein position of genomic position
        protein_position = position_in_cds // 3 + 1 if position_in_cds % 3 != 0 else position_in_cds // 3
        # Get positions of codon --> in case the mutations starts after the codon
        codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
        delins = False
        first_aa = ""
        last_aa = ""
        last_position = ""
        codon_last = []
        ref_aa = ""

        if position_in_cds % 3 == 0:
            protein_position += 1
            codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
            first_aa = str(Seq(self.genome[
                              (protein_start - 1) + codon[0] - 1: (protein_start - 1) + codon[2]]
                              ).translate())
            # If the deletion is longer than one codon get last codon in deletion
            if deletion_length / 3 > 1 and deletion_length % 3 == 0:
                last_position = protein_position + deletion_length // 3 - 1
                codon_last = [last_position * 3 - 2, last_position * 3 - 1, last_position * 3]
                last_aa = str(Seq(self.genome[
                                 (protein_start - 1) + codon_last[0] - 1: (protein_start - 1) + codon_last[2]]
                                 ).translate())

        # Deletion starts at first base of codon (includes second and third)
        elif position_in_cds % 3 == 1:
            # Created deleted and WT sequence and check if AA has changed -->
            # If the codon has not changed the following AA are affected
            wt_start = (protein_start - 1) + codon[0]
            wt_end = (protein_start - 1) + codon[2]
            wt_aa = str(Seq(self.genome[wt_start - 1: wt_end]).translate())

            mt_end = wt_start + deletion_length
            mt_aa = str(Seq(self.genome[wt_start - 1] + self.genome[mt_end: mt_end + 2]).translate())
            # If affected codon doesn't change the following codons are affected
            if wt_aa == mt_aa:
                protein_position += 1
                codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
                first_aa = str(Seq(self.genome[
                                  (protein_start - 1) + codon[0] - 1: (protein_start - 1) + codon[2]]
                                  ).translate())
                # Get last codon and AA of deletion
                if deletion_length / 3 > 1 and deletion_length % 3 == 0:
                    last_position = protein_position + deletion_length // 3 - 1
                    codon_last = [last_position * 3 - 2, last_position * 3 - 1, last_position * 3]
                    last_aa = str(Seq(self.genome[
                                     (protein_start - 1) + codon_last[0] - 1: (protein_start - 1) + codon_last[2]]
                                     ).translate())
            else:
                first_aa = wt_aa
                # Get last position and AA for delins case
                if deletion_length / 3 >= 1 and deletion_length % 3 == 0:
                    delins = True
                    inserted_aa = mt_aa
                    last_position = protein_position + deletion_length // 3
                    codon_last = [last_position * 3 - 2, last_position * 3 - 1, last_position * 3]
                    last_aa = str(Seq(self.genome[
                                     (protein_start - 1) + codon_last[0] - 1: (protein_start - 1) + codon_last[2]]
                                     ).translate())

        # Deletion starts at second base of codon (includes third base)
        elif position_in_cds % 3 == 2:
            wt_start = (protein_start - 1) + codon[0]
            wt_end = (protein_start - 1) + codon[2]
            wt_aa = str(Seq(self.genome[wt_start - 1: wt_end]).translate())

            mt_end = (protein_start - 1) + codon[1] + deletion_length
            mt_aa = str(Seq(self.genome[wt_start - 1] + self.genome[wt_start + 1] + self.genome[mt_end]).translate())
            ref_aa = str(Seq(self.genome[wt_start - 1: mt_end + 1]).translate())
            if mt_aa == wt_aa:
                # The next AA is affected
                protein_position += 1
                codon = [protein_position * 3 - 2, protein_position * 3 - 1, protein_position * 3]
                first_aa = str(Seq(self.genome[
                                   (protein_start - 1) + codon[0] - 1: (protein_start - 1) + codon[2]]
                                   ).translate())
                # Get last codon and AA of deletion
                if deletion_length / 3 > 1 and deletion_length % 3 == 0:
                    last_position = protein_position + deletion_length // 3 - 1
                    codon_last = [last_position * 3 - 2, last_position * 3 - 1, last_position * 3]
                    last_aa = str(Seq(self.genome[
                                     (protein_start - 1) + codon_last[0] - 1: (protein_start - 1) + codon_last[2]]
                                     ).translate())
            else:
                first_aa = wt_aa
                # Get last position and AA for delins case
                if deletion_length / 3 >= 1 and deletion_length % 3 == 0:
                    delins = True
                    inserted_aa = mt_aa
                    last_position = protein_position + deletion_length // 3
                    codon_last = [last_position * 3 - 2, last_position * 3 - 1, last_position * 3]
                    last_aa = str(Seq(self.genome[
                                      (protein_start - 1) + codon_last[0] - 1: (protein_start - 1) + codon_last[2]]
                                      ).translate())

        hgvs = ""
        if deletion_length % 3 != 0:
            hgvs = "p.{aa}{pos}fs".format(
                aa=first_aa,
                pos=protein_position
            )
            ref_aa = first_aa
        elif delins:
            hgvs = "p.{first_aa}{first_pos}_{last_aa}{last_pos}delins{ins}".format(
                first_aa=first_aa,
                first_pos=protein_position,
                last_aa=last_aa,
                last_pos=last_position,
                ins=inserted_aa
            )
            ref_aa = str(Seq(self.genome[(protein_start - 1) + codon[0] - 1:
                                         (protein_start - 1) + codon_last[2]]).translate())
        else:
            if deletion_length / 3 > 1 and deletion_length % 3 == 0:
                hgvs = "p.{first_aa}{first_pos}_{last_aa}{last_pos}del".format(
                    first_aa=first_aa,
                    first_pos=protein_position,
                    last_aa=last_aa,
                    last_pos=last_position
                )
                ref_aa = str(Seq(self.genome[(protein_start - 1) + codon[0] - 1:
                                             (protein_start - 1) + codon_last[2]]).translate())
            else:
                hgvs = "p.{aa}{pos}del".format(
                    aa=first_aa,
                    pos=protein_position
                )
                ref_aa = first_aa
        return hgvs, ref_aa, protein_position

    def translate_snps_into_protein(self, variant, protein_coordinates):
        position = variant.get("position")
        protein_coordinates = protein_coordinates.query("start <= @genomic_position & end >= @genomic_position").get(
            ["gene", "start", "end"])
        protein_start = protein_coordinates.get("start").item()


    def _parse_mutation_sites(self, variant_string: str):
        """
        Parse mutation list from pango constellation files. This code
        is based on the parser found in the utility script of scorpio
        for now limited to SNPs, Indels and Insertion in protein coordinates.

        s:Y144- --> Deletion on AA level
        del:21765:6 --> Deletion on nucleotide level

        nuc:C16176T  --> SNV on nucleotide level
        s:A570D --> SNV on AA level
        n:RG203KR --> MNV on AA level

        nuc:22205+GAGCCAGAA --> Insertion on nucleotide level

        """
        snp_pattern = re.compile('([ACTGUN]+)([0-9]+)([ACTGUN]+)')
        insertion_pattern = re.compile(r'(\w+):(\d+)\+([a-zA-Z]+)')
        aa_mutation = re.compile(r'([a-zA-Z-*]+)(\d+)([a-zA-Z-*]*)')
        mutation_list = []
        # This parsing is simplified code based on the utility script of Scorpio
        this_mut = variant_string.split(":")
        # Insertions -> eiter nucleotide or protein coordinates
        if "+" in variant_string:
            match = re.match(insertion_pattern, variant_string)
            if not match:
                logger.warning("Could not parse the following site definition: {}".format(variant_string))
                return [None]
            # Does not work with aa level insertions but None are present in the current constellation files
            mutation_info = {"site": variant_string,
                             "type": "Insertion",
                             "position": match[2],
                             "length": len(match[2]),
                             "protein": self._match_protein_name(match[1]),
                             "alternate": match[3],
                             "level": "nuc",
                             "hgvs_p": None}
            if not this_mut[0] in ["snp", "nuc"]:
                mutation_info["level"] = "aa"
            mutation_list.append(mutation_info)
        # Deletions -> coordinates on nucleotide level
        elif this_mut[0] == "del":
            length = int(this_mut[2])
            protein = self.find_gene(this_mut[1])
            mutation_info = {"site": variant_string,
                             "variant_id": "",
                             "type": "DELETION",
                             "genomic_position": this_mut[1],
                             "position": "",
                             "reference": "",
                             "alternate": "del",
                             "length": length,
                             "protein": protein,
                             "level": "aa",
                             "hgvs_p": None}
            hgvs, ref, pos_aa = self.get_hgvs_from_nuc_deletion(mutation_info)
            mutation_info["hgvs_p"] = hgvs
            mutation_info["reference"] = ref
            mutation_info["position"] = pos_aa
            mutation_info["variant_id"] = "{}:{}{}{}".format(mutation_info["protein"], mutation_info["reference"],
                                                             mutation_info["position"], mutation_info["alternate"])
            mutation_list.append(mutation_info)
        # SNPs -> coordinates on nucleotide level, non-synonymous or intergenic?
        elif this_mut[0] in ["snp", "nuc"]:
            match = re.match(snp_pattern, this_mut[1])
            if not match:
                logger.warning("Could not parse the following site definition: {}".format(variant_string))
                return [None]
            mutation_info = {"site": variant_string,
                             "type": "SNV",
                             "position": int(match[2]),
                             "reference": match[1],
                             "alternate": match[3],
                             "protein": None,
                             "level": "nuc",
                             "hgvs_p": None}
            mutation_info['variant_id'] = "{}:{}{}{}".format(mutation_info["protein"], mutation_info["reference"],
                                                             mutation_info["position"], mutation_info["alternate"])
            mutation_list.append(mutation_info)
        # Amino acid substitutions --> can be SNVs, MNVs, Deletions and Indels
        else:
            match = re.match(aa_mutation, this_mut[1])
            if not match:
                logger.warning("Could not parse the following site definition: {}".format(variant_string))
                return [None]
            protein = self._match_protein_name(this_mut[0])
            reference = match[1]
            alternate = match[3]
            position = int(match[2])
            if alternate not in ['-', 'del']:
                if len(reference) == len(alternate) and len(reference) > 1:
                    for protein_pos, aa_mutations in enumerate(zip(reference, alternate), start=position):
                        variant_id = "{}:{}{}{}".format(protein, aa_mutations[0], protein_pos, aa_mutations[1])
                        mutation_info = {"site": variant_id,
                                         "variant_id": variant_id,
                                         "hgvs_p": "p.{}{}{}".format(aa_mutations[0], protein_pos, aa_mutations[1]),
                                         "type": "SNV",
                                         "position": protein_pos,
                                         "reference": aa_mutations[0],
                                         "alternate": aa_mutations[1],
                                         "protein": protein,
                                         "ambiguous_alternate": False,
                                         "length": len(aa_mutations[0]),
                                         "level": "aa"}
                        mutation_list.append(mutation_info)
                else:
                    variant_id = "{}:{}{}{}".format(protein, reference, position, alternate)
                    mutation_info = {"site": variant_string,
                                     "variant_id": variant_id,
                                     "hgvs_p": "p.{}{}{}".format(reference, position, alternate),
                                     "type": "SNV",
                                     "position": position,
                                     "reference": reference,
                                     "alternate": alternate,
                                     "protein": protein,
                                     "ambiguous_alternate": True if alternate == "" else False,
                                     "length": len(reference),
                                     "level": "aa"}
                    mutation_list.append(mutation_info)
            else:
                alternate = "del" if alternate == "-" else "del"
                variant_id = "{}:{}{}{}".format(protein, reference, position, alternate)
                if len(reference) > 1:
                    hgvs_p = "p.{}{}_{}{}del".format(reference[0], position, reference[-1],
                                                     position + len(reference)-1)
                else:
                    hgvs_p = "p.{}{}del".format(reference, position)

                mutation_info = {"site": variant_string,
                                 "variant_id": variant_id,
                                 "hgvs_p": hgvs_p,
                                 "type": "DELETION",
                                 "position": position,
                                 "reference": reference,
                                 "alternate": alternate,
                                 "protein": protein,
                                 "ambiguous_alternate": False,
                                 "length": len(reference),
                                 "level": "aa"}
                mutation_list.append(mutation_info)
        return mutation_list

    @staticmethod
    def _create_constellation_pango_mapping(lineage_constellation: dict):
        """
        Creates a dictionary, mapping pango lineage identifiers to constellation labels. This mapping
        is required in find_parent_sites to find the constellation label for a parent lineage id as
        parent lineages are given as pango ids.
        """
        lineage_mapping = {}
        for constellation, annotation in lineage_constellation.items():
            for pango_id in annotation["pangolin_lineage_list"]:
                lineage_mapping[pango_id] = constellation
        return lineage_mapping

    def _find_parent_sites(self, lineage: str, lineage_constellation: dict, pango_constellation_mapping: dict):
        """
        Include for each constellation also the mutations present up in the phylogenetic tree.
        """
        # Returns parent as pangolin identifier
        parent_pango = lineage_constellation[lineage].get("parent_lineage_id", None)
        # Get corresponding constellation label to look up mutations
        parent = pango_constellation_mapping.get(parent_pango, None)
        sites = lineage_constellation[lineage].get("lineage_mutations", None)
        if parent is None:
            return sites
        else:
            sites.extend([x for x in self._find_parent_sites(parent, lineage_constellation, pango_constellation_mapping)
                          if x not in sites])
        return sites

    def _find_constellation_files(self):
        """
        Glob all lineage constellation files found in constellation repository
        """
        json_files = list(Path(self.lineage_constellation_directory).glob('*.json'))
        return json_files

    def _read_constellation_files(self):
        """
        Parse lineage defining constellation JSON files and store information
        """
        lineage_constellation = {}
        # Iterate over all json constellation files
        for this_file in self._find_constellation_files():
            with open(this_file, "r") as fd:
                data = json.load(fd)
            # Get lineage annotation information
            pangolin_lineage_list = set(data["variant"].get("Pango_lineages", []))
            if not pangolin_lineage_list:
                # MRCA lineage info is either not present or literally an empty string ""
                lineage_name = data["variant"].get("mrca_lineage", "")
                if not lineage_name:
                    lineage_name = data["variant"].get("lineage_name", "")
                pangolin_lineage_list = set([lineage_name])
            constellation_label = data["label"]
            # Optional records not present in all lineage constellations
            who_label = data["variant"].get("WHO_label", None)
            phe_label = data["variant"].get("PHE_label", None)
            voc_date = None
            vui_date = None
            variant_of_concern = False
            variant_under_investigation = False
            tags = "|".join(data["tags"])
            # Extract VOC/VOI/VUI information from tags field
            for this_tag in data["tags"]:
                if this_tag.startswith("VOC-") or this_tag.startswith("VOC "):
                    variant_of_concern = True
                    if this_tag.startswith("VOC-"):
                        voc_date = datetime.strptime(this_tag.lstrip("VOC-"), "%y%b-%d")
                        voc_date = voc_date.strftime("%Y-%m-%d")
                    # edge case for B1.1.7 constellation file
                    else:
                        voc_date = datetime.strptime(this_tag.lstrip("VOC "), "%Y%m/%d")
                        voc_date = voc_date.strftime("%Y-%m-%d")
                # Some constellations bring both VUI and V tags, others only one of the two.
                # Therefore, we have to check both cases, even if this means that we overwrite the information once.
                if this_tag.startswith("V-"):
                    variant_under_investigation = True
                    vui_date = datetime.strptime(this_tag.lstrip("V-"), "%y%b-%d")
                    vui_date = vui_date.strftime("%Y-%m-%d")
                if this_tag.startswith("VUI-"):
                    variant_under_investigation = True
                    vui_date = datetime.strptime(this_tag.lstrip("VUI-"), "%y%b-%d")
                    vui_date = vui_date.strftime("%Y-%m-%d")

            parent_lineage_id = data["variant"].get("parent_lineage", None)
            # Drop incompatible lineage calls from pangolin identifier list
            # These calls are incompatible with the scorpio constellation
            incompatible_lineage_calls = set(data["variant"].get("incompatible_lineage_calls", set()))
            pangolin_lineage_list = pangolin_lineage_list - incompatible_lineage_calls
            lineage_mutations = [self._parse_mutation_sites(x) for x in data["sites"]]
            lineage_mutations = [mut for x in lineage_mutations for mut in x if mut is not None]
            lineage_constellation[constellation_label] = {
                "pangolin_lineage_list": pangolin_lineage_list,
                "who_label": who_label,
                "phe_label": phe_label,
                "voc_date": voc_date,
                "vui_date": vui_date,
                "variant_of_concern": variant_of_concern,
                "variant_under_investigation": variant_under_investigation,
                "tags": tags,
                "parent_lineage_id": parent_lineage_id,
                "lineage_mutations": lineage_mutations
            }
        return lineage_constellation

    def _fill_lineage_mutation_table(self, lineage_constellation: dict):
        """
        Store constellation information in database
        """
        count_lineages = 0
        count_mutations = 0
        all_mutations = []
        seen_mutations = set()

        # Get a list of pangolin lineage ids with a PHE label
        lineages_with_phe = [y.get("pangolin_lineage_list") for x, y in lineage_constellation.items()
                             if y.get("phe_label")]
        lineages_with_phe = {x for lineage in lineages_with_phe for x in lineage}
        for constellation, annotation in lineage_constellation.items():
            # Drop "duplicated" lineages, as we do not need to annotate them twice in the database.
            # This is the case when a separate constellation file was created for an already existing lineage,
            # just to store it with an additional mutation, e.g. B.1.1.7+E484K.json and B.1.617.2+K417N.
            # Because it is still the same lineage, we store in the table only the constellation with a given PHE label
            pango_lineages = annotation["pangolin_lineage_list"]
            if not annotation["phe_label"]:
                pango_lineages = annotation["pangolin_lineage_list"] - lineages_with_phe

            for pangolin_lineage_id in pango_lineages:
                # Multiple pangolin identifiers can be presented in a constellation.
                # That is the case, e.g., for the delta sub-lineages AY.1 and AY.2, that were created to track local
                # spread. For now, we store the pango-lineages together with the constellation label as unique keys
                lineage = Lineages(
                    pango_lineage_id=pangolin_lineage_id,
                    constellation_id=constellation,
                    who_label=annotation["who_label"],
                    phe_label=annotation["phe_label"],
                    parent_lineage_id=annotation["parent_lineage_id"],
                    voc_date=annotation["voc_date"],
                    vui_date=annotation["vui_date"],
                    variant_under_investigation=annotation["variant_under_investigation"],
                    variant_of_concern=annotation["variant_of_concern"],
                    tags=annotation["tags"]
                )
                self.session.add(lineage)
                self.session.commit()
                count_lineages += 1
            # Store mutations in list of seen mutations
            for this_mut in annotation["lineage_mutations"]:
                if this_mut["level"] != "aa":
                    continue
                if this_mut["variant_id"] not in seen_mutations:
                    seen_mutations.add(this_mut["variant_id"])
                    all_mutations.append(this_mut)
        logger.info("Loaded into the database {} lineages".format(count_lineages))

        # Store mutations in database
        for this_mut in all_mutations:
            # Skip variants on nucleotide level for now
            mutation = LineageDefiningVariants(
                variant_id=this_mut["variant_id"],
                hgvs=this_mut["hgvs_p"],
                variant_type=this_mut["type"],
                protein=this_mut["protein"],
                position=this_mut["position"],
                reference=this_mut["reference"],
                alternate=this_mut["alternate"],
                ambiguous_alternate=this_mut["ambiguous_alternate"]
            )
            self.session.add(mutation)
            self.session.commit()
            count_mutations += 1
        logger.info("Loaded into the database {} mutations".format(count_mutations))

    def _fill_relation_ship_table(self, lineage_constellation):
        """
        Store lineage/mutation N:M relationships
        """
        count_relationship = 0
        mutation_lineage_mapping = {}
        pango_constellation_mapping = self._create_constellation_pango_mapping(lineage_constellation)
        # Pick up all mutations/lineage combinations
        for constellation, annotation in lineage_constellation.items():
            # Include mutations from parent lineages as well
            all_mutations_of_lineage = self._find_parent_sites(constellation, lineage_constellation,
                                                               pango_constellation_mapping)
            for this_mut in all_mutations_of_lineage:
                # Skip mutations that are not on AA level
                # This can be safely deleted as soon as nucleotide level problems are solved
                if this_mut["level"] != "aa":
                    continue
                var_id = this_mut["variant_id"]
                if var_id not in mutation_lineage_mapping.keys():
                    mutation_lineage_mapping[var_id] = annotation["pangolin_lineage_list"]
                else:
                    # Create union of pango lineages with the mutation
                    mutation_lineage_mapping[var_id] = mutation_lineage_mapping[var_id] | annotation["pangolin_lineage_list"]

        for variant_id, lineages in mutation_lineage_mapping.items():
            for pango_id in lineages:
                lineage_variant = LineageVariant(
                    pango_lineage_id=pango_id,
                    variant_id=variant_id
                )
                self.session.add(lineage_variant)
                self.session.commit()
                count_relationship += 1
        logger.info("Loaded into the database {} lineage-variant relationships".format(count_relationship))

    def load_data(self):
        # Fill lineage annotation table
        lineage_constellation = self._read_constellation_files()
        self._fill_lineage_mutation_table(lineage_constellation)
        self._fill_relation_ship_table(lineage_constellation)
