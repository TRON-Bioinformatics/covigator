import os
import json
from pathlib import Path

import pandas as pd
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

    def load_data(self):
        count_lineages = 0
        # Iterate over all json constellation files
        for this_file in self._find_constellation_files():
            with open(this_file, "r") as fd:
                data = json.load(fd)
            # Get lineage annotation information
            pangolin_lineage_list = set(data["variant"]["Pango_lineages"])
            constellation_label = data["label"]
            # Optional records not present in all lineage constellations
            who_label = data["variant"].get("WHO_label", "")
            phe_label = data["variant"].get("PHE_label", "")
            date = ""
            variant_of_concern = False
            variant_under_investigation = False
            tags = data["tags"]
            for this_tag in tags:
                if this_tag.startswith("VOC-"):
                    variant_of_concern = True
                    # Remove prefix from phe label date
                    constellation_date = phe_label.lstrip("VOC-")
                elif this_tag.startswith("V-"):
                    variant_under_investigation = True
                    constellation_date = phe_label.lstrip("V-")
            #if phe_label:
            #    # The PHE label shows the date when the variant was classified as variant of concern/interest
            #    # Check if string begins with VOC- or just V- and removes prefix to get classification date
            #    if phe_label.startswith("VOC-"):
            #        variant_of_concern = True
            #        # Remove prefix from phe label date
            #        constellation_date = phe_label.lstrip("VOC-")
            #    else:
            #        variant_of_interest = True
            #        constellation_date = phe_label.lstrip("V-")

                # Dates are formatted in the json file the following way: 2021APR-02
            constellation_date = datetime.strptime(constellation_date, "%y%b-%d")
            date = constellation_date.strftime("%Y-%m-%d")

            parent_lineage_id = data["variant"].get("parent_lineage", "")
            incompatible_lineage_calls = set(data["variant"].get("incompatible_lineage_calls", []))
            # Drop incompatible lineage calls from pango identifier list
            pangolin_lineage_list = pangolin_lineage_list - incompatible_lineage_calls
            for pangolin_lineage_id in pangolin_lineage_list:
                # It seems we can have multiple pangolin identifiers in the pango designation field.
                # That is the case for delta sublineages that were classified AY.X
                # I decide for now to store pango lineages together with the constellation label as unique keys

                lineage = Lineages(
                    pangolin_lineage_id=pangolin_lineage_id,
                    constellation_label=constellation_label,
                    who_label=who_label,
                    phe_label=phe_label,
                    parent_lineage_id=parent_lineage_id,
                    date=date,
                    variant_under_investigation=variant_under_investigation,
                    variant_of_concern=variant_of_concern
                )
                self.session.add(lineage)
                self.session.commit()
                count_lineages += 1
        logger.info("Loaded into the database {} lineages".format(count_lineages))
