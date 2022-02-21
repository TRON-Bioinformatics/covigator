import pandas as pd
from logzero import logger
from sqlalchemy import func, desc
from sqlalchemy.orm import Session
from covigator import SYNONYMOUS_VARIANT
from covigator.database.model import DataSource, PrecomputedOccurrence
from covigator.database.queries import Queries


NUMBER_TOP_OCCURRENCES = 1000


class TopOccurrencesLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)

    def load(self):

        database_rows = self._read_top_occurrent_variants()

        # delete all rows before starting
        self.session.query(PrecomputedOccurrence).delete()
        self.session.commit()

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedOccurrence.__tablename__))

    def _read_top_occurrent_variants(self):
        # gets the top occurrent variants for each source
        top_occurring_variants_ena = None
        try:
            top_occurring_variants_ena = self.get_top_occurring_variants(
                top=NUMBER_TOP_OCCURRENCES, source=DataSource.ENA.name)
        except ValueError as e:
            logger.exception(e)
            logger.error("No top occurrences for ENA data")
        top_occurring_variants_gisaid = None
        try:
            top_occurring_variants_gisaid = self.get_top_occurring_variants(
                top=NUMBER_TOP_OCCURRENCES, source=DataSource.GISAID.name)
        except ValueError:
            logger.error("No top occurrences for GISAID data")
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
        return database_rows

    def _row_to_top_occurrence(self, row, source):
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

    def get_top_occurring_variants(self, top, source: str):
        klass = self.queries.get_variant_observation_klass(source)
        query = self.session.query(
            klass.variant_id, klass.hgvs_p, klass.gene_name, klass.pfam_name,
            klass.annotation_highest_impact, func.count().label('total')) \
            .filter(klass.annotation_highest_impact != SYNONYMOUS_VARIANT)
        query = query.group_by(klass.variant_id, klass.hgvs_p, klass.gene_name,
                      klass.pfam_name, klass.annotation_highest_impact) \
            .order_by(desc('total')).limit(top)
        top_occurring_variants = pd.read_sql(query.statement, self.session.bind)

        # calculate frequency
        count_samples = self.queries.count_samples(source=source)
        top_occurring_variants['frequency'] = top_occurring_variants.total.transform(
            lambda x: round(float(x) / count_samples, 3))

        # add counts for every month
        top_occurring_variants = self._get_counts_per_month(
            top_occurring_variants=top_occurring_variants, source=source)

        return top_occurring_variants

    def _get_counts_per_month(self, top_occurring_variants, source: str):
        variant_counts_by_month = []
        for _, variant in top_occurring_variants.iterrows():
            variant_counts_by_month.append(self.queries.get_variant_counts_by_month(variant.variant_id, source=source))
        if len(variant_counts_by_month) > 1:
            top_occurring_variants_by_month = pd.concat(variant_counts_by_month)
            # get total count of samples per month to calculate the frequency by month
            sample_counts_by_month = self.queries.get_sample_counts_by_month(source=source)
            top_occurring_variants_by_month = pd.merge(
                left=top_occurring_variants_by_month, right=sample_counts_by_month, how="left", on="month")
            top_occurring_variants_by_month["frequency_by_month"] = \
                (top_occurring_variants_by_month["count"] / top_occurring_variants_by_month["sample_count"]). \
                    transform(lambda x: round(x, 3))
            # join both tables with total counts and counts per month
            top_occurring_variants = pd.merge(
                left=top_occurring_variants, right=top_occurring_variants_by_month, on="variant_id", how="left")
            # format the month column appropriately
            top_occurring_variants.month = top_occurring_variants.month.transform(
                lambda d: "{}-{:02d}".format(d.year, int(d.month)))
        return top_occurring_variants
