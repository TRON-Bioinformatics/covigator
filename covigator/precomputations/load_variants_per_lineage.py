import pandas as pd
from logzero import logger
from sqlalchemy import func, and_
from sqlalchemy.orm import Session
from covigator import SYNONYMOUS_VARIANT
from covigator.database.model import DataSource, PrecomputedVariantsPerLineage
from covigator.database.queries import Queries


class VariantsPerLineageLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)

    def load(self):

        database_rows = self._read_variants_per_lineage()

        # delete all rows before starting
        self.session.query(PrecomputedVariantsPerLineage).delete()
        self.session.commit()

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedVariantsPerLineage.__tablename__))

    def _read_variants_per_lineage(self):
        # gets the top occurrent variants for each source
        variants_per_lineage_ena = None
        try:
            variants_per_lineage_ena = self.get_variants_per_lineage(source=DataSource.ENA.name)
        except ValueError as e:
            logger.exception(e)
            logger.error("No top occurrences for ENA data")
        variants_per_lineage_gisaid = None
        try:
            variants_per_lineage_gisaid = self.get_variants_per_lineage(source=DataSource.GISAID.name)
        except ValueError:
            logger.error("No top occurrences for GISAID data")
        database_rows = []
        # stores the precomputed data
        if variants_per_lineage_ena is not None:
            for index, row in variants_per_lineage_ena.iterrows():
                # add entries per gene
                database_rows.append(self._row_to_variants_per_lineage(row, source=DataSource.ENA))
        if variants_per_lineage_gisaid is not None:
            for index, row in variants_per_lineage_gisaid.iterrows():
                # add entries per gene
                database_rows.append(self._row_to_variants_per_lineage(row, source=DataSource.GISAID))
        return database_rows

    def _row_to_variants_per_lineage(self, row, source):
        return PrecomputedVariantsPerLineage(
            lineage=row["lineage"],
            variant_id=row["variant_id"],
            country=row["country"],
            count_observations=row["count_observations"],
            source=source
        )

    def get_variants_per_lineage(self, source: str):
        klass_variant_observation = self.queries.get_variant_observation_klass(source)
        klass_sample = self.queries.get_sample_klass(source)
        query = self.session.query(
            klass_variant_observation.variant_id,
            klass_sample.country,
            klass_sample.pangolin_lineage.label("lineage"),
            func.count().label('count_observations'))\
            .filter(and_(
                klass_variant_observation.annotation_highest_impact != SYNONYMOUS_VARIANT,
                klass_sample.pangolin_lineage != None,
                klass_sample.pangolin_lineage != "")) \
            .join(klass_sample, klass_variant_observation.sample == klass_sample.run_accession) \
            .group_by(klass_sample.pangolin_lineage,
                      klass_sample.country,
                      klass_variant_observation.variant_id)
        variants_per_lineage = pd.read_sql(query.statement, self.session.bind)
        return variants_per_lineage
