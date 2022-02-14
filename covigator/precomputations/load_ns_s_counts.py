from typing import List

import pandas as pd
from logzero import logger
from sqlalchemy.orm import Session
from covigator import SYNONYMOUS_VARIANT, MISSENSE_VARIANT
from covigator.database.model import DataSource, PrecomputedSynonymousNonSynonymousCounts, RegionType
from covigator.database.queries import Queries


class NsSCountsLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)

    def load(self):

        counts_by_gene = self._load_counts(region=RegionType.GENE)
        counts_by_domain = self._load_counts(region=RegionType.DOMAIN)
        database_rows = counts_by_gene + counts_by_domain

        # delete all rows before starting
        self.session.query(PrecomputedSynonymousNonSynonymousCounts).delete()
        self.session.commit()

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows),
                                                    PrecomputedSynonymousNonSynonymousCounts.__tablename__))

    def _load_counts(self, region: RegionType) -> List[PrecomputedSynonymousNonSynonymousCounts]:
        data_s_ena = self._count_variant_observations_by_source_annotation_and_region(
            source=DataSource.ENA, annotation=SYNONYMOUS_VARIANT, region=region)
        data_ns_ena = self._count_variant_observations_by_source_annotation_and_region(
            source=DataSource.ENA, annotation=MISSENSE_VARIANT, region=region)
        data_s_gisaid = self._count_variant_observations_by_source_annotation_and_region(
            source=DataSource.GISAID, annotation=SYNONYMOUS_VARIANT, region=region)
        data_ns_gisaid = self._count_variant_observations_by_source_annotation_and_region(
            source=DataSource.GISAID, annotation=MISSENSE_VARIANT, region=region)

        data_ena = pd.merge(
            left=data_s_ena, right=data_ns_ena, on=["month", "region_name", "country"], how='outer').fillna(0)
        data_gisaid = pd.merge(
            left=data_s_gisaid, right=data_ns_gisaid, on=["month", "region_name", "country"], how='outer').fillna(0)

        database_rows = []
        if data_ena is not None:
            database_rows.extend(self._dataframe_to_model(
                data=data_ena, source=DataSource.ENA, region=region))
        if data_gisaid is not None:
            database_rows.extend(self._dataframe_to_model(
                data=data_gisaid, source=DataSource.GISAID, region=region))

        return database_rows

    def _dataframe_to_model(self, data, source: DataSource, region: RegionType):

        database_rows = []
        for index, row in data.iterrows():
            database_rows.append(self._row_to_model(row, source, region))

        if region == RegionType.GENE:
            # if we are processing genes then we calculate also for the coding region overall
            coding_region_ns = {}
            coding_region_s = {}
            for index, row in data.iterrows():
                month = row['month']
                country = row['country']
                ns = row["ns"]
                s = row["s"]
                coding_region_ns[(month, country)] = coding_region_ns.get((month, country), 0) + ns
                coding_region_s[(month, country)] = coding_region_s.get((month, country), 0) + s
            for month, country in coding_region_ns.keys():
                database_rows.append(PrecomputedSynonymousNonSynonymousCounts(
                    month=month,
                    region_type=RegionType.CODING_REGION,
                    country=country,
                    ns=coding_region_ns.get((month, country), 0),
                    s=coding_region_s.get((month, country), 0),
                    source=source
                ))

        return database_rows

    def _row_to_model(self, row, source, region: RegionType):
        # add entries per gene
        return PrecomputedSynonymousNonSynonymousCounts(
            month=row['month'],
            region_type=region,
            region_name=row["region_name"],
            country=row['country'],
            ns=row["ns"],
            s=row["s"],
            source=source
        )

    def _count_variant_observations_by_source_annotation_and_region(
            self, source: DataSource, annotation: str, region: RegionType):

        variant_observation_table = self.queries.get_variant_observation_klass(source=source.name).__tablename__
        sample_table = self.queries.get_sample_klass(source=source.name).__tablename__

        sql_query = """
                select count(*) as {count_name}, date_trunc('month', s.collection_date) as month, 
                    vo.{region_field} as region_name, s.country 
                from {variant_observation_table} as vo join {sample_table} as s on vo.sample = s.run_accession 
                where vo.annotation_highest_impact = '{annotation}'
                    and s.collection_date is not null and vo.{region_field} is not null
                group by date_trunc('month', s.collection_date), vo.{region_field}, s.country;
                """.format(variant_observation_table=variant_observation_table,
                           sample_table=sample_table,
                           region_field="gene_name" if region == RegionType.GENE else "pfam_name",
                           count_name="s" if annotation == SYNONYMOUS_VARIANT else "ns",
                           annotation=annotation)
        data = pd.read_sql_query(sql_query, self.session.bind)
        return data
