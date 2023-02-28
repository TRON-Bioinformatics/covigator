import os
import json
import re
from pathlib import Path

from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session
from covigator.database.model import Lineages
from logzero import logger
from datetime import datetime


class LineageAnnotationsLoader:

    LINEAGE_CONSTELLATION_DIRECTORY = "constellations/constellations/definitions"

    def __init__(self, session: Session):
        self.session = session

        self.lineage_constellation_directory = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.LINEAGE_CONSTELLATION_DIRECTORY)

    def _find_constellation_files(self):
        # Glob all lineage constellation files from subdirectory
        json_files = list(Path(self.lineage_constellation_directory).glob('*.json'))
        return json_files

    def _match_protein_name(self, protein_name):
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

    def _parse_mutation_sites(self, site_list):
        """
        Parse mutation list from pango constellation file
        for now limited to SNPs in genomic coordinates
        """
        snp_pattern = re.compile('([ACTGUN]+)([0-9]+)([ACTGUN]+)')
        insertion_pattern = re.compile(r'(\w+):(\d+)\+([a-zA-Z]+)')
        aa_mutation = re.compile(r'([a-zA-Z-*]+)(\d+)([a-zA-Z-*]*)')
        mutation_list = []
        for this_site in site_list:
            # This parsing is simplified code based on the utility script of Scorpio
            this_mut = this_site.split(":")
            # Insertions -> eiter nucleotide or protein coordinates
            if "+" in this_mut[1]:
                match = re.match(insertion_pattern, this_mut[1])
                if not match:
                    logger.warning("Could not parse the following site definition: {}".format(this_site))
                    continue
                # ToDo support insertions and ambigous mutations that are in nuc coordinates
                # How to get reference allele so we can compare?
                mutation_info = {"site": this_site, "type": "Insertion", "position": match[2], "length": length(match[2]),
                                 "protein": self._match_protein_name(match[1]), "alternate": match[3]}
                # Not typed in scorpio --> how to handle?2
            # Deletions -> coordinates in nucleotide coordinates
            elif this_mut[0] == "del":
                length = int(this_mut[2])
                mutation_info = {"site": this_site, "type": "Deletion", "position": this_mut[1], "length": length,
                                 "protein": None}
            # SNPs -> coordinates in nucleotide coordinates
            elif this_mut[0] in ["snp", "nuc"]:
                match = re.match(snp_pattern, this_mut[1])
                if not match:
                    logger.warning("Could not parse the following site definition: {}".format(this_site))
                    continue
                mutation_info = {"site": this_site, "type": "SNV", "position": int(match[2]),
                                 "reference": match[1], "alternate": match[3], "protein": None}
            # Amino acid substitutions --> can be SNVs, Deletions and Indels?
            # Ambigous codes are problematic for the DB
            else:
                match = re.match(aa_mutation, this_mut[1])
                if not match:
                    logger.warning("Could not parse the following site definition: {}".format(this_site))
                    continue
                protein = self._match_protein_name(this_mut[0])
                mutation_info = {"site": this_site, "type": "", "position": int(match[2]),
                                 "reference": match[1], "alternate": match[3], "protein": protein
                                 }
                if mutation_info["alternate"] == '':
                    mutation_info["ambiguous_alternate"] = True
                elif mutation_info["alternate"] in ["-", "del"]:
                    mutation_info["type"] = "AA_Deletion"
                    mutation_info["length"] = len(mutation_info["reference"])
                    mutation_info["ambiguous_alternate"] = False
                else:
                    mutation_info["ambiguous_alternate"] = False
            mutation_list.append(mutation_info)

        return mutation_list

    def _read_constellation_files(self):
        lineage_constellation = {}
        # Iterate over all json constellation files
        for this_file in self._find_constellation_files():
            with open(this_file, "r") as fd:
                data = json.load(fd)
            # Get lineage annotation information
            pangolin_lineage_list = set(data["variant"].get("Pango_lineages", ""))
            if not pangolin_lineage_list:
                pangolin_lineage_list = set(data["variant"].get("mrca_lineage"))
            constellation_label = data["label"]
            # Optional records not present in all lineage constellations
            who_label = data["variant"].get("WHO_label", "")
            phe_label = data["variant"].get("PHE_label", "")
            date = None
            variant_of_concern = False
            variant_under_investigation = False
            tags = "|".join(data["tags"])
            sites = self._parse_mutation_sites(data['sites'])
            if phe_label:
                if phe_label.startswith("VOC-"):
                    variant_of_concern = True
                    constellation_date = phe_label.lstrip("VOC-")
                if phe_label.startswith("V-"):
                    variant_under_investigation = True
                    constellation_date = phe_label.lstrip("V-")

                constellation_date = datetime.strptime(constellation_date, "%y%b-%d")
                date = constellation_date.strftime("%Y-%m-%d")
            # Dates are formatted in the json file the following way: 2021APR-02

            parent_lineage_id = data["variant"].get("parent_lineage", "")
            incompatible_lineage_calls = set(data["variant"].get("incompatible_lineage_calls", []))
            # Drop incompatible lineage calls from pango identifier list
            pangolin_lineage_list = pangolin_lineage_list - incompatible_lineage_calls
            lineage_constellation[constellation_label] = {
                "pangolin_lineage_list": pangolin_lineage_list,
                "who_label": who_label,
                "phe_label": phe_label,
                "date": date,
                "variant_of_concern": variant_of_concern,
                "variant_under_investigation": variant_under_investigation,
                "tags": tags,
                "parent_lineage_id": parent_lineage_id,
                "sites": sites
            }
        return lineage_constellation

    def _find_parent_sites(self, lineage, lineage_constellation):
        """
        Include for each constellation also the mutations present up in the phylogenetic tree.
        """
        parent = lineage_constellation[lineage].get("parent_lineage", "")
        sites = lineage_constellation[lineage].get("sites")
        # Convert sites into set for quick union operation
        if not parent:
            return sites
        else:
            sites.extend([x for x in self._find_parent_sites(parent, lineage_constellation) if x not in sites])
        return sites

    def _fill_lineage_table(self, lineage_constellation):
        count_lineages = 0
        lineages_with_phe = [y.get("pangolin_lineage_list") for x, y in lineage_constellation.items() if y.get("phe_label")]
        for constellation, annotation in lineage_constellation.items():
            # Drop "duplicated" lineages as we don't need to annotate them twice in the database.
            # This is the case when a separate constellation file was created for an already existing lineage
            # just to store it with an additional mutation. As it is still the same lineage we only store in the lineage
            # table the constellation with a set PHE label
            if not annotation["phe_label"] and annotation["pangolin_lineage_list"] in lineages_with_phe:
                continue
            for pangolin_lineage_id in annotation["pangolin_lineage_list"]:
                # It seems we can have multiple pangolin identifiers in the pango designation field.
                # That is the case for delta sublineages AY.1 adn AY.2
                # I decide for now to store pango lineages together with the constellation label as unique keys
                lineage = Lineages(
                    pango_lineage_id=pangolin_lineage_id,
                    constellation_id=constellation,
                    who_label=annotation["who_label"],
                    phe_label=annotation["phe_label"],
                    parent_lineage_id=annotation["parent_lineage_id"],
                    date=annotation["date"],
                    variant_under_investigation=annotation["variant_under_investigation"],
                    variant_of_concern=annotation["variant_of_concern"],
                    tags=annotation["tags"]
                )
                self.session.add(lineage)
                self.session.commit()
                count_lineages += 1
        logger.info("Loaded into the database {} lineages".format(count_lineages))

    def _fill_sites_table(self, lineage_constellation):
        pass


    def load_data(self):
        # Fill lineage annotation table
        lineage_constellation = self._read_constellation_files()
        self._fill_lineage_table(lineage_constellation)
