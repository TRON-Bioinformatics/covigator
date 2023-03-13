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
        """
        Glob all lineage constellation files found in constellation repository
        """
        json_files = list(Path(self.lineage_constellation_directory).glob('*.json'))
        return json_files

    def _read_constellation_files(self):
        """
        Parse lineage defining constellation files and store information
        as a dict of dicts
        """
        lineage_constellation = {}
        # Iterate over all json constellation files
        for this_file in self._find_constellation_files():
            with open(this_file, "r") as fd:
                data = json.load(fd)
            # Get lineage annotation information
            pangolin_lineage_list = set(data["variant"].get("Pango_lineages", []))
            if not pangolin_lineage_list:
                # MRCA lineages are either not present or an empty string ""
                lineage_name = data["variant"].get("mrca_lineage", "")
                if not lineage_name:
                    lineage_name = data["variant"].get("lineage_name", "")
                pangolin_lineage_list = set([lineage_name])
            constellation_label = data["label"]
            # Optional records not present in all lineage constellations
            who_label = data["variant"].get("WHO_label", "")
            phe_label = data["variant"].get("PHE_label", "")
            voc_date = None
            vui_date = None
            variant_of_concern = False
            variant_under_investigation = False
            tags = "|".join(data["tags"])
            for this_tag in data["tags"]:
                if this_tag.startswith("VOC-"):
                    variant_of_concern = True
                    voc_date = datetime.strptime(this_tag.lstrip("VOC-"), "%y%b-%d")
                    voc_date = voc_date.strftime("%Y-%m-%d")
                if this_tag.startswith("V-"):
                    variant_under_investigation = True
                    vui_date = datetime.strptime(this_tag.lstrip("V-"), "%y%b-%d")
                    vui_date = vui_date.strftime("%Y-%m-%d")
            parent_lineage_id = data["variant"].get("parent_lineage", "")
            # Drop incompatible lineage calls from pangolin identifier list
            # These calls are used by
            incompatible_lineage_calls = set(data["variant"].get("incompatible_lineage_calls", []))
            pangolin_lineage_list = pangolin_lineage_list - incompatible_lineage_calls
            lineage_constellation[constellation_label] = {
                "pangolin_lineage_list": pangolin_lineage_list,
                "who_label": who_label,
                "phe_label": phe_label,
                "voc_date": voc_date,
                "vui_date": vui_date,
                "variant_of_concern": variant_of_concern,
                "variant_under_investigation": variant_under_investigation,
                "tags": tags,
                "parent_lineage_id": parent_lineage_id
            }
        return lineage_constellation

    def _fill_lineage_table(self, lineage_constellation):
        """
        Store constellation information in database
        """
        count_lineages = 0
        # Get a list of pangolin lineage ids with a PHE label
        lineages_with_phe = [y.get("pangolin_lineage_list") for x, y in lineage_constellation.items() if y.get("phe_label")]
        lineages_with_phe = {x for lineage in lineages_with_phe for x in lineage}
        for constellation, annotation in lineage_constellation.items():
            # Drop "duplicated" lineages as we don't need to annotate them twice in the database.
            # This is the case when a separate constellation file was created for an already existing lineage
            # just to store it with an additional mutation. As it is still the same lineage we only store in the lineage
            # table the constellation with a set PHE label
            pango_lineages = annotation["pangolin_lineage_list"]
            if not annotation["phe_label"]:
                # Skip lineages that occur multiple times in the constellation files
                pango_lineages = annotation["pangolin_lineage_list"] - lineages_with_phe

            for pangolin_lineage_id in pango_lineages:
                # Multiple pangolin identifiers can point to the same lineage.
                # That is the case for e.g. the delta sub-lineages AY.1 and AY.2 that were created to track local
                # spreadings
                # For now we store the pango lineages together with the constellation label as unique keys
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
        logger.info("Loaded into the database {} lineages".format(count_lineages))

    def load_data(self):
        # Fill lineage annotation table
        lineage_constellation = self._read_constellation_files()
        self._fill_lineage_table(lineage_constellation)
