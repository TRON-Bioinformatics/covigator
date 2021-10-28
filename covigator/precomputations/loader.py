import pandas as pd
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.dashboard.tabs.recurrent_mutations import BIN_SIZE_VALUES
from covigator.database.database import Database
from covigator.database.model import PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, \
    VARIANT_OBSERVATION_TABLE_NAME, PrecomputedIndelLength, \
    PrecomputedAnnotation, VariantObservation, DataSource, \
    PrecomputedTableCounts, Variant, \
    SubclonalVariantObservation, Sample, PrecomputedVariantAbundanceHistogram
from logzero import logger

from covigator.database.queries import Queries
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.precomputations.load_top_occurrences import TopOccurrencesLoader


class PrecomputationsLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)
        self.ns_s_counts_loader = NsSCountsLoader(session=session)
        self.top_occurrences_loader = TopOccurrencesLoader(session=session)

    def load(self):
        logger.info("Starting precomputations...")
        self.load_table_counts()
        logger.info("Done with table counts (1/8)")
        self.load_variant_abundance_histogram()
        logger.info("Done with variant abundance histogram (2/8)")
        self.load_counts_variants_per_sample()
        logger.info("Done with count variants per sample (3/8)")
        self.load_count_substitutions()
        logger.info("Done with count base subsitutions (4/8)")
        self.load_indel_length()
        logger.info("Done with indel length (5/8)")
        self.load_annotation()
        logger.info("Done with effect annotations (6/8)")
        self.top_occurrences_loader.load()
        logger.info("Done with top occurrent variants (7/8)")
        self.ns_s_counts_loader.load()
        logger.info("Done with NS S counts (8/8)")

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
            if histogram is not None:
                histogram["bin_size"] = bin_size
                histogram["source"] = None
                histograms.append(histogram)

        for bin_size in BIN_SIZE_VALUES:
            histogram = self.queries.get_variant_abundance_histogram(
                bin_size=bin_size, source=DataSource.ENA.name, cache=False)
            if histogram is not None:
                histogram["bin_size"] = bin_size
                histogram["source"] = DataSource.ENA
                histograms.append(histogram)

        for bin_size in BIN_SIZE_VALUES:
            histogram = self.queries.get_variant_abundance_histogram(
                bin_size=bin_size, source=DataSource.GISAID.name, cache=False)
            if histogram is not None:
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
    precomputer = PrecomputationsLoader(session=database.get_database_session())
    #precomputer.load_dn_ds()
    #precomputer.load_table_counts()
    #precomputer.load_variant_abundance_histogram()
    precomputer.load_top_occurrences()
