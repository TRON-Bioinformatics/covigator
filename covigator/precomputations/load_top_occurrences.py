from typing import List

import pandas as pd
from logzero import logger
from sqlalchemy.orm import Session
from covigator import SYNONYMOUS_VARIANT, MISSENSE_VARIANT
from covigator.database.model import DataSource, PrecomputedSynonymousNonSynonymousCounts, RegionType, \
    VARIANT_OBSERVATION_TABLE_NAME, SAMPLE_ENA_TABLE_NAME, SAMPLE_GISAID_TABLE_NAME, PrecomputedOccurrence
from covigator.database.queries import Queries


NUMBER_TOP_OCCURRENCES = 1000


class TopOccurrencesLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)

    def load(self):

        # gets the top occurrent variants for each source and overall
        top_occurring_variants_ena = None
        try:
            top_occurring_variants_ena = self.queries.get_top_occurring_variants(
                top=NUMBER_TOP_OCCURRENCES, source=DataSource.ENA.name)
        except ValueError as e:
            logger.exception(e)
            logger.error("No top occurrences for ENA data")

        top_occurring_variants_gisaid = None
        try:
            top_occurring_variants_gisaid = self.queries.get_top_occurring_variants(
                top=NUMBER_TOP_OCCURRENCES, source=DataSource.GISAID.name)
        except ValueError:
            logger.error("No top occurrences for GISAID data")

        top_occurring_variants = None
        try:
            top_occurring_variants = self.queries.get_top_occurring_variants(top=NUMBER_TOP_OCCURRENCES)
        except ValueError:
            logger.error("No top occurrences")

        # delete all rows before starting
        self.session.query(PrecomputedOccurrence).delete()
        self.session.commit()

        database_rows = []
        # stores the precomputed data
        if top_occurring_variants_ena is not None:
            for index, row in top_occurring_variants_ena.iterrows():
                # add entries per gene
                database_rows.append(self._row_to_top_occurrence(row, source=DataSource.ENA))

        if top_occurring_variants_gisaid is not None:
            for index, row in top_occurring_variants_gisaid.iterrows():
                # add entries per gene
                database_rows.append(self._row_to_top_occurrence(row, source=DataSource.GISAID))

        if top_occurring_variants is not None:
            for index, row in top_occurring_variants.iterrows():
                # add entries per gene
                database_rows.append(self._row_to_top_occurrence(row))

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedOccurrence.__tablename__))

    def _row_to_top_occurrence(self, row, source=None):
        return PrecomputedOccurrence(
            total=row["total"],
            frequency=row["frequency"],
            variant_id=row["variant_id"],
            hgvs_p=row["hgvs_p"],
            gene_name=row["gene_name"],
            domain=row["pfam_name"],
            annotation=row["annotation_highest_impact"],
            source=source,
            month=row["month"],
            count=row["count"],
            frequency_by_month=row["frequency_by_month"],
        )
