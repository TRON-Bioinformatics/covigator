import os

import sqlalchemy
from sqlalchemy.orm import Session
from covigator.database.model import Conservation
from logzero import logger
import pandas as pd


class ConservationLoader:

    CONSERVATION_FILENAME = "wuhCor1.mutDepletionConsHMM.bed"
    CONSERVATION_SARBECOVIRUS_FILENAME = "wuhCor1.mutDepletionSarbecovirusConsHMM.bed"
    CONSERVATION_VERTEBRATES_FILENAME = "wuhCor1.mutDepletionVertebrateCoVConsHMM.bed"

    def __init__(self, session: Session):

        # these files were downloaded from https://github.com/ernstlab/ConsHMM_CoV
        self.conservation_filename = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.CONSERVATION_FILENAME)
        self.conservation_sarbecovirus_filename = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.CONSERVATION_SARBECOVIRUS_FILENAME)
        self.conservation_vertebrates_filename = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.CONSERVATION_VERTEBRATES_FILENAME)

        self.session = session

    def load_data(self):
        # read conservation data from BED files
        conservation = pd.read_csv(
            self.conservation_filename, skiprows=1, names=["chromosome", "start", "end", "conservation"], sep="\t")
        conservation_sarbecovirus = pd.read_csv(
            self.conservation_sarbecovirus_filename, skiprows=1,
            names=["chromosome", "start", "end", "conservation_sarbecovirus"], sep="\t")
        conservation_vertebrates = pd.read_csv(
            self.conservation_vertebrates_filename, skiprows=1,
            names=["chromosome", "start", "end", "conservation_vertebrates"], sep="\t")

        # joins conservation data in a single table
        data = conservation.set_index(["chromosome", "start", "end"]).join(
            conservation_sarbecovirus.set_index(["chromosome", "start", "end"])
        ).join(
            conservation_vertebrates.set_index(["chromosome", "start", "end"])
        )
        data.reset_index(inplace=True)
        data.fillna(0, inplace=True)

        # persist the conservation table
        metadata = sqlalchemy.schema.MetaData(bind=self.session.bind, reflect=True)
        conservation_table = sqlalchemy.Table(Conservation.__tablename__, metadata, autoload=True)
        self.session.execute(conservation_table.insert(), data.to_dict(orient='records'))
        self.session.commit()

        logger.info("Loaded into the database {} conservation data points".format(data.shape[0]))
