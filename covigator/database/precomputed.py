import pandas as pd
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.dashboard.tabs.variants import BIN_SIZE_VALUES
from covigator.database.database import Database
from covigator.database.model import PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, \
    VARIANT_OBSERVATION_TABLE_NAME, PrecomputedIndelLength, \
    PrecomputedAnnotation, VariantObservation, DataSource, PrecomputedOccurrence, PrecomputedDnDs, \
    SAMPLE_ENA_TABLE_NAME, SAMPLE_GISAID_TABLE_NAME, PrecomputedDnDsByDomain, PrecomputedTableCounts, Variant, \
    SubclonalVariantObservation, Sample, PrecomputedVariantAbundanceHistogram
from logzero import logger

from covigator.database.queries import Queries

NUMBER_TOP_OCCURRENCES = 1000


class Precomputer:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)

    def load_counts_variants_per_sample(self):

        # reads the data from the database
        sql_query = """
        select count(*) as count, counts.number_mutations, counts.source, counts.variant_type, counts.gene_name from (
            select count(*) as number_mutations, source, variant_type, gene_name, sample from {table_name}
            group by source, variant_type, gene_name, sample)
            as counts group by counts.number_mutations, counts.source, counts.variant_type, counts.gene_name;
            """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, counts.number_mutations, counts.source, counts.variant_type from (
            select count(*) as number_mutations, source, variant_type, sample from {table_name}
            group by source, variant_type, sample)
            as counts group by counts.number_mutations, counts.source, counts.variant_type;
            """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        # delete all rows before starting
        self.session.query(PrecomputedVariantsPerSample).delete()
        self.session.commit()

        # stores the precomputed data
        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedVariantsPerSample(
                count=row["count"],
                number_mutations=row["number_mutations"],
                source=row["source"],
                variant_type=row["variant_type"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))

        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedVariantsPerSample(
                count=row["count"],
                number_mutations=row["number_mutations"],
                source=row["source"],
                variant_type=row["variant_type"]
            ))

        self.session.add_all(database_rows)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedVariantsPerSample.__tablename__))

    def load_count_substitutions(self):

        sql_query = """
        select * from (select count(*) as count, reference, alternate, source, variant_type, gene_name 
            from {table_name}
            group by reference, alternate, source, variant_type, gene_name
            order by count desc) 
            as counts where counts.count > 10;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select * from (select count(*) as count, reference, alternate, source, variant_type 
            from {table_name}
            group by reference, alternate, source, variant_type
            order by count desc) 
            as counts where counts.count > 10;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        # delete all rows before starting
        self.session.query(PrecomputedSubstitutionsCounts).delete()
        self.session.commit()

        # stores the precomputed data
        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedSubstitutionsCounts(
                count=row["count"],
                reference=row["reference"],
                alternate=row["alternate"],
                source=row["source"],
                variant_type=row["variant_type"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))

        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedSubstitutionsCounts(
                count=row["count"],
                reference=row["reference"],
                alternate=row["alternate"],
                source=row["source"],
                variant_type=row["variant_type"]
            ))

        self.session.add_all(database_rows)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedSubstitutionsCounts.__tablename__))

    def load_indel_length(self):

        sql_query = """
        select count(*) as count, length, source, gene_name 
            from {table_name}
            where variant_type != 'SNV'
            group by length, source, gene_name
            order by count asc;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, length, source, gene_name 
            from {table_name}
            where variant_type != 'SNV'
            group by length, source, gene_name
            order by count desc;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        # delete all rows before starting
        self.session.query(PrecomputedIndelLength).delete()
        self.session.commit()

        # stores the precomputed data
        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedIndelLength(
                count=row["count"],
                length=row["length"],
                source=row["source"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))

        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedIndelLength(
                count=row["count"],
                length=row["length"],
                source=row["source"]
            ))

        self.session.add_all(database_rows)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedIndelLength.__tablename__))

    def load_annotation(self):

        sql_query = """
        select count(*) as count, annotation_highest_impact, source, gene_name 
            from {table_name}
            where annotation_highest_impact is not null
            group by annotation_highest_impact, source, gene_name
            order by count asc;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, annotation_highest_impact, source 
            from {table_name}
            where annotation_highest_impact is not null
            group by annotation_highest_impact, source
            order by count asc;
        """.format(table_name=VARIANT_OBSERVATION_TABLE_NAME)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        # delete all rows before starting
        self.session.query(PrecomputedAnnotation).delete()
        self.session.commit()

        # stores the precomputed data
        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedAnnotation(
                count=row["count"],
                annotation=row["annotation_highest_impact"],
                source=row["source"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))

        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedAnnotation(
                count=row["count"],
                annotation=row["annotation_highest_impact"],
                source=row["source"]
            ))

        self.session.add_all(database_rows)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedAnnotation.__tablename__))

    def load_top_occurrences(self):

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
            annotation=row["annotation_highest_impact"],
            source=source,
            month=row["month"],
            count=row["count"],
            frequency_by_month=row["frequency_by_month"],
        )

    def load_dn_ds(self):
        sql_query_ds_ena = """
                select count(*) as ds, date_trunc('month', s.collection_date) as month, vo.gene_name, s.country 
                from {variant_observation_table} as vo join {sample_ena_table} as s on vo.sample = s.run_accession 
                where vo.source = 'ENA' and vo.annotation_highest_impact = 'synonymous_variant' 
                    and s.collection_date is not null
                group by date_trunc('month', s.collection_date), vo.gene_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_ena_table=SAMPLE_ENA_TABLE_NAME)
        data_ds_ena = pd.read_sql_query(sql_query_ds_ena, self.session.bind)

        sql_query_dn_ena = """
                select count(*) as dn, date_trunc('month', s.collection_date) as month, vo.gene_name, s.country 
                from {variant_observation_table} as vo join {sample_ena_table} as s on vo.sample = s.run_accession 
                where vo.source = 'ENA' and vo.annotation_highest_impact = 'missense_variant' 
                    and s.collection_date is not null
                group by date_trunc('month', s.collection_date), vo.gene_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_ena_table=SAMPLE_ENA_TABLE_NAME)
        data_dn_ena = pd.read_sql_query(sql_query_dn_ena, self.session.bind)

        data_ena = pd.merge(left=data_ds_ena, right=data_dn_ena, on=["month", "gene_name", "country"], how='outer').fillna(0)

        sql_query_ds_gisaid = """
                select count(*) as ds, date_trunc('month', s.date) as month, vo.gene_name, s.country 
                from {variant_observation_table} as vo join {sample_gisaid_table} as s on vo.sample = s.run_accession 
                    and s.date is not null
                where vo.source = 'GISAID' and vo.annotation_highest_impact = 'synonymous_variant' 
                group by date_trunc('month', s.date), vo.gene_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_gisaid_table=SAMPLE_GISAID_TABLE_NAME)
        data_ds_gisaid = pd.read_sql_query(sql_query_ds_gisaid, self.session.bind)

        sql_query_dn_gisaid = """
                select count(*) as dn, date_trunc('month', s.date) as month, vo.gene_name, s.country 
                from {variant_observation_table} as vo join {sample_gisaid_table} as s on vo.sample = s.run_accession 
                where vo.source = 'GISAID' and vo.annotation_highest_impact = 'missense_variant'
                    and s.date is not null 
                group by date_trunc('month', s.date), vo.gene_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_gisaid_table=SAMPLE_GISAID_TABLE_NAME)
        data_dn_gisaid = pd.read_sql_query(sql_query_dn_gisaid, self.session.bind)

        data_gisaid = pd.merge(left=data_ds_gisaid, right=data_dn_gisaid, on=["month", "gene_name", "country"],
                            how='outer').fillna(0)

        # delete all rows before starting
        self.session.query(PrecomputedDnDs).delete()
        self.session.commit()

        database_rows = []
        if data_ena is not None:
            for index, row in data_ena.iterrows():
                # add entries per gene
                database_rows.append(PrecomputedDnDs(
                    month=row["month"],
                    gene_name=row["gene_name"],
                    country=row["country"],
                    dn=row["dn"],
                    ds=row["ds"],
                    source=DataSource.ENA
                ))
        if data_gisaid is not None:
            for index, row in data_gisaid.iterrows():
                # add entries per gene
                database_rows.append(PrecomputedDnDs(
                    month=row["month"],
                    gene_name=row["gene_name"],
                    country=row["country"],
                    dn=row["dn"],
                    ds=row["ds"],
                    source=DataSource.GISAID
                ))

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedDnDs.__tablename__))


    def load_dn_ds_by_domain(self):
        sql_query_ds_ena = """
                select count(*) as ds, date_trunc('month', s.collection_date) as month, vo.pfam_name, s.country 
                from {variant_observation_table} as vo join {sample_ena_table} as s on vo.sample = s.run_accession 
                where vo.source = 'ENA' and vo.annotation_highest_impact = 'synonymous_variant' 
                    and vo.pfam_name is not null and s.collection_date is not null
                group by date_trunc('month', s.collection_date), vo.pfam_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_ena_table=SAMPLE_ENA_TABLE_NAME)
        data_ds_ena = pd.read_sql_query(sql_query_ds_ena, self.session.bind)

        sql_query_dn_ena = """
                select count(*) as dn, date_trunc('month', s.collection_date) as month, vo.pfam_name, s.country 
                from {variant_observation_table} as vo join {sample_ena_table} as s on vo.sample = s.run_accession 
                where vo.source = 'ENA' and vo.annotation_highest_impact = 'missense_variant' 
                    and vo.pfam_name is not null and s.collection_date is not null
                group by date_trunc('month', s.collection_date), vo.pfam_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_ena_table=SAMPLE_ENA_TABLE_NAME)
        data_dn_ena = pd.read_sql_query(sql_query_dn_ena, self.session.bind)

        data_ena = pd.merge(left=data_ds_ena, right=data_dn_ena, on=["month", "pfam_name", "country"], how='outer').fillna(0)

        sql_query_ds_gisaid = """
                select count(*) as ds, date_trunc('month', s.date) as month, vo.pfam_name, s.country 
                from {variant_observation_table} as vo join {sample_gisaid_table} as s on vo.sample = s.run_accession 
                where vo.source = 'GISAID' and vo.annotation_highest_impact = 'synonymous_variant' 
                    and vo.pfam_name is not null and s.date is not null
                group by date_trunc('month', s.date), vo.pfam_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_gisaid_table=SAMPLE_GISAID_TABLE_NAME)
        data_ds_gisaid = pd.read_sql_query(sql_query_ds_gisaid, self.session.bind)

        sql_query_dn_gisaid = """
                select count(*) as dn, date_trunc('month', s.date) as month, vo.pfam_name, s.country 
                from {variant_observation_table} as vo join {sample_gisaid_table} as s on vo.sample = s.run_accession 
                where vo.source = 'GISAID' and vo.annotation_highest_impact = 'missense_variant' 
                    and vo.pfam_name is not null and s.date is not null
                group by date_trunc('month', s.date), vo.pfam_name, s.country;
                """.format(variant_observation_table=VARIANT_OBSERVATION_TABLE_NAME,
                           sample_gisaid_table=SAMPLE_GISAID_TABLE_NAME)
        data_dn_gisaid = pd.read_sql_query(sql_query_dn_gisaid, self.session.bind)

        data_gisaid = pd.merge(left=data_ds_gisaid, right=data_dn_gisaid, on=["month", "pfam_name", "country"],
                            how='outer').fillna(0)

        # delete all rows before starting
        self.session.query(PrecomputedDnDs).delete()
        self.session.commit()

        database_rows = []
        if data_ena is not None:
            for index, row in data_ena.iterrows():
                # add entries per gene
                database_rows.append(PrecomputedDnDsByDomain(
                    month=row["month"],
                    domain=row["pfam_name"],
                    country=row["country"],
                    dn=row["dn"],
                    ds=row["ds"],
                    source=DataSource.ENA
                ))
        if data_gisaid is not None:
            for index, row in data_gisaid.iterrows():
                # add entries per gene
                database_rows.append(PrecomputedDnDsByDomain(
                    month=row["month"],
                    domain=row["pfam_name"],
                    country=row["country"],
                    dn=row["dn"],
                    ds=row["ds"],
                    source=DataSource.GISAID
                ))

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedDnDsByDomain.__tablename__))

    def load_table_counts(self):

        count_variants = self.queries.count_variants(cache=False)
        count_samples = self.queries.count_samples(cache=False)
        count_samples_ena = self.queries.count_samples(source=DataSource.ENA.name, cache=False)
        count_samples_gisaid = self.queries.count_samples(source=DataSource.GISAID.name, cache=False)
        count_variant_observations = self.queries.count_variant_observations(cache=False)
        count_variant_observations_ena = self.queries.count_variant_observations(source=DataSource.ENA.name, cache=False)
        count_variant_observations_gisaid = self.queries.count_variant_observations(source=DataSource.GISAID.name, cache=False)
        count_subclonal_variant_observations = self.queries.count_subclonal_variant_observations(cache=False)
        count_subclonal_variant_unique = self.queries.count_unique_subclonal_variant(cache=False)
        count_subclonal_variant_unique_only_subclonal = self.queries.count_unique_only_subclonal_variant(cache=False)
        count_countries = self.queries.count_countries(cache=False)
        count_countries_ena = self.queries.count_countries(source=DataSource.ENA.name, cache=False)
        count_countries_gisaid = self.queries.count_countries(source=DataSource.GISAID.name, cache=False)

        # delete all rows before starting
        self.session.query(PrecomputedTableCounts).delete()
        self.session.commit()

        database_rows = [
            PrecomputedTableCounts(table=Variant.__name__, count=count_variants),
            PrecomputedTableCounts(table=VariantObservation.__name__, count=count_variant_observations),
            PrecomputedTableCounts(
                table=VariantObservation.__name__, count=count_variant_observations_ena,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(
                table=VariantObservation.__name__, count=count_variant_observations_gisaid,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__, count=count_subclonal_variant_observations),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__ + "_unique", count=count_subclonal_variant_unique),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__ + "_unique_only_subclonal",
                count=count_subclonal_variant_unique_only_subclonal),
            PrecomputedTableCounts(table=Sample.__name__, count=count_samples),
            PrecomputedTableCounts(table=Sample.__name__, count=count_samples_ena,
                                   factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(table=Sample.__name__, count=count_samples_gisaid,
                                   factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
            PrecomputedTableCounts(table=PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY, count=count_countries),
            PrecomputedTableCounts(table=PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY, count=count_countries_ena,
                                   factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(table=PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY, count=count_countries_gisaid,
                                   factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
        ]

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedTableCounts.__tablename__))

    def load_variant_abundance_histogram(self):

        histograms = []
        for bin_size in BIN_SIZE_VALUES:
            histogram = self.queries.get_variant_abundance_histogram(bin_size=bin_size, cache=False)
            histogram["bin_size"] = bin_size
            histogram["source"] = None
            histograms.append(histogram)

        for bin_size in BIN_SIZE_VALUES:
            histogram = self.queries.get_variant_abundance_histogram(
                bin_size=bin_size, source=DataSource.ENA.name, cache=False)
            histogram["bin_size"] = bin_size
            histogram["source"] = DataSource.ENA
            histograms.append(histogram)

        for bin_size in BIN_SIZE_VALUES:
            histogram = self.queries.get_variant_abundance_histogram(
                bin_size=bin_size, source=DataSource.GISAID.name, cache=False)
            histogram["bin_size"] = bin_size
            histogram["source"] = DataSource.GISAID
            histograms.append(histogram)

        # delete all rows before starting
        self.session.query(PrecomputedVariantAbundanceHistogram).delete()
        self.session.commit()

        database_rows = []
        if len(histograms) > 0:
            for index, row in pd.concat(histograms).iterrows():
                # add entries per gene
                database_rows.append(PrecomputedVariantAbundanceHistogram(
                    position_bin=row["position_bin"],
                    count_unique_variants=row["count_unique_variants"],
                    count_variant_observations=row["count_variant_observations"],
                    bin_size=row["bin_size"],
                    source=row["source"],
                ))

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(
            len(database_rows), PrecomputedVariantAbundanceHistogram.__tablename__))


if __name__ == '__main__':
    database = Database(initialize=True, config=Configuration())
    precomputer = Precomputer(session=database.get_database_session())
    #precomputer.load_dn_ds_by_domain()
    #precomputer.load_dn_ds()
    #precomputer.load_table_counts()
    #precomputer.load_variant_abundance_histogram()
    precomputer.load_top_occurrences()
