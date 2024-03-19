import os
import pandas as pd
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session
from covigator.database.model import IntrahostPrediction
from logzero import logger


class IntrahostPredictionLoader:

    PREDICTIONS = "all_mutations_predictions.parquet"

    def __init__(self, session: Session):
        self.session = session

        self.intrahost_predictions = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.PREDICTIONS)

    def load_data(self):
        # reads the parquet file into dataframe
        data = pd.read_parquet(self.intrahost_predictions)

        data = self._remove_duplicated_genes(data)
        
        count_mutations = 0
        for mut in data.itertuples():
            try:
                IntrahostPrediction(
                        variant_id = mut.variant_id,
                        gene_name = mut.gene_name,
                        referene = mut.reference_amino_acid,
                        alternate = mut.alternate_amino_acid,
                        clonal_dnds = mut.clonal_dnds,
                        subclonal_dnds = mut.subclonal_dnds,
                        sift_score = mut.sift_score,
                        prediction = mut.prediction,
                        prediction_probability = mut.prediction_proba
                )
                self.session.add(gene)
                self.session.commit()
            except IntegrityError:
                logger.warn("Mutation {} is duplicated in input".format(mut.variant_id))               
                self.session.rollback()
            count_mutations += 1
        logger.info("Loaded into the database {} intrahost model predictions".format(count_mutations))
