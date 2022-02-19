import pandas as pd
from sqlalchemy.orm import Session
from covigator.dashboard.tabs.recurrent_mutations import BIN_SIZE_VALUES
from covigator.database.model import PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, \
    PrecomputedIndelLength, PrecomputedAnnotation, VariantObservation, DataSource, \
    PrecomputedTableCounts, Variant, SubclonalVariantObservation, PrecomputedVariantAbundanceHistogram, \
    SampleEna, SampleGisaid, GisaidVariantObservation, GisaidVariant
from logzero import logger

from covigator.database.queries import Queries
from covigator.precomputations.load_ns_s_counts import NsSCountsLoader
from covigator.precomputations.load_top_occurrences import TopOccurrencesLoader
from covigator.precomputations.load_variants_per_lineage import VariantsPerLineageLoader


class PrecomputationsLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)
        self.ns_s_counts_loader = NsSCountsLoader(session=session)
        self.top_occurrences_loader = TopOccurrencesLoader(session=session)
        self.variants_per_lineage_loader = VariantsPerLineageLoader(session=session)

    def load(self):
        logger.info("Starting precomputations...")
        self.load_table_counts()
        logger.info("Done with table counts (1/9)")
        self.load_variant_abundance_histogram()
        logger.info("Done with variant abundance histogram (2/9)")
        self.load_counts_variants_per_sample()
        logger.info("Done with count variants per sample (3/9)")
        self.load_count_substitutions()
        logger.info("Done with count base subsitutions (4/9)")
        self.load_indel_length()
        logger.info("Done with indel length (5/9)")
        self.load_annotation()
        logger.info("Done with effect annotations (6/9)")
        self.top_occurrences_loader.load()
        logger.info("Done with top occurrent variants (7/9)")
        self.ns_s_counts_loader.load()
        logger.info("Done with NS S counts (8/9)")
        self.variants_per_lineage_loader.load()
        logger.info("Done with variants per lineage (9/9)")

    def load_counts_variants_per_sample(self):

        # reads the data from the database
        database_rows_ena = self._read_count_variants_per_sample(data_source=DataSource.ENA)
        database_rows_gisaid = self._read_count_variants_per_sample(data_source=DataSource.GISAID)

        # delete all rows before starting
        self.session.query(PrecomputedVariantsPerSample).delete()
        self.session.commit()

        self.session.add_all(database_rows_ena + database_rows_gisaid)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows_ena) + len(database_rows_gisaid),
                                                    PrecomputedVariantsPerSample.__tablename__))

    def _read_count_variants_per_sample(self, data_source: DataSource):

        table_name = self.queries.get_variant_observation_klass(source=data_source.name).__tablename__

        sql_query = """
        select count(*) as count, counts.number_mutations, counts.variant_type, counts.gene_name from (
            select count(*) as number_mutations, variant_type, gene_name, sample from {table_name}
            group by variant_type, gene_name, sample)
            as counts group by counts.number_mutations, counts.variant_type, counts.gene_name;
            """.format(table_name=table_name)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, counts.number_mutations, counts.variant_type from (
            select count(*) as number_mutations, variant_type, sample from {table_name}
            group by variant_type, sample)
            as counts group by counts.number_mutations, counts.variant_type;
            """.format(table_name=table_name)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        database_rows = []
        for index, row in data.iterrows():
            database_rows.append(PrecomputedVariantsPerSample(
                count=row["count"],
                number_mutations=row["number_mutations"],
                source=data_source.name,
                variant_type=row["variant_type"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))
        for index, row in data_without_gene.iterrows():
            database_rows.append(PrecomputedVariantsPerSample(
                count=row["count"],
                number_mutations=row["number_mutations"],
                source=data_source.name,
                variant_type=row["variant_type"]
            ))
        return database_rows

    def load_count_substitutions(self):

        database_rows_ena = self._read_count_substitutions(data_source=DataSource.ENA)
        database_rows_gisaid = self._read_count_substitutions(data_source=DataSource.GISAID)

        # delete all rows before starting
        self.session.query(PrecomputedSubstitutionsCounts).delete()
        self.session.commit()

        self.session.add_all(database_rows_ena + database_rows_gisaid)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows_ena) + len(database_rows_gisaid),
                                                    PrecomputedSubstitutionsCounts.__tablename__))

    def _read_count_substitutions(self, data_source: DataSource):

        table_name = self.queries.get_variant_observation_klass(source=data_source.name).__tablename__

        sql_query = """
        select * from (select count(*) as count, reference, alternate, variant_type, gene_name 
            from {table_name}
            group by reference, alternate, variant_type, gene_name
            order by count desc) 
            as counts where counts.count > 10;
        """.format(table_name=table_name)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select * from (select count(*) as count, reference, alternate, variant_type 
            from {table_name}
            group by reference, alternate, variant_type
            order by count desc) 
            as counts where counts.count > 10;
        """.format(table_name=table_name)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        database_rows = []
        for index, row in data.iterrows():
            database_rows.append(PrecomputedSubstitutionsCounts(
                count=row["count"],
                reference=row["reference"],
                alternate=row["alternate"],
                source=data_source.name,
                variant_type=row["variant_type"],
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))
        for index, row in data_without_gene.iterrows():
            database_rows.append(PrecomputedSubstitutionsCounts(
                count=row["count"],
                reference=row["reference"],
                alternate=row["alternate"],
                source=data_source.name,
                variant_type=row["variant_type"]
            ))
        return database_rows

    def load_indel_length(self):

        database_rows_ena = self._read_indel_length(data_source=DataSource.ENA)
        database_rows_gisaid = self._read_indel_length(data_source=DataSource.GISAID)

        # delete all rows before starting
        self.session.query(PrecomputedIndelLength).delete()
        self.session.commit()

        self.session.add_all(database_rows_ena + database_rows_gisaid)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows_ena) + len(database_rows_gisaid),
                                                    PrecomputedIndelLength.__tablename__))

    def _read_indel_length(self, data_source: DataSource):

        table_name = self.queries.get_variant_observation_klass(source=data_source.name).__tablename__

        sql_query = """
        select count(*) as count, length, gene_name 
            from {table_name}
            where variant_type != 'SNV'
            group by length, gene_name
            order by count asc;
        """.format(table_name=table_name)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, length, gene_name 
            from {table_name}
            where variant_type != 'SNV'
            group by length, gene_name
            order by count desc;
        """.format(table_name=table_name)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)

        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedIndelLength(
                count=row["count"],
                length=row["length"],
                source=data_source.name,
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))
        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedIndelLength(
                count=row["count"],
                length=row["length"],
                source=data_source.name
            ))
        return database_rows

    def load_annotation(self):

        database_rows_ena = self._read_annotations(data_source=DataSource.ENA)
        database_rows_gisaid = self._read_annotations(data_source=DataSource.GISAID)

        # delete all rows before starting
        self.session.query(PrecomputedAnnotation).delete()
        self.session.commit()

        self.session.add_all(database_rows_ena + database_rows_gisaid)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows_ena) + len(database_rows_gisaid), PrecomputedAnnotation.__tablename__))

    def _read_annotations(self, data_source: DataSource):

        table_name = self.queries.get_variant_observation_klass(source=data_source.name).__tablename__

        sql_query = """
        select count(*) as count, annotation_highest_impact, gene_name 
            from {table_name}
            where annotation_highest_impact is not null
            group by annotation_highest_impact, gene_name
            order by count asc;
        """.format(table_name=table_name)
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*) as count, annotation_highest_impact 
            from {table_name}
            where annotation_highest_impact is not null
            group by annotation_highest_impact
            order by count asc;
        """.format(table_name=table_name)
        data_without_gene = pd.read_sql_query(sql_query, self.session.bind)
        # stores the precomputed data
        database_rows = []
        for index, row in data.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedAnnotation(
                count=row["count"],
                annotation=row["annotation_highest_impact"],
                source=data_source.name,
                gene_name=row["gene_name"] if row["gene_name"] is not None else "intergenic"
            ))
        for index, row in data_without_gene.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedAnnotation(
                count=row["count"],
                annotation=row["annotation_highest_impact"],
                source=data_source.name
            ))
        return database_rows

    def load_table_counts(self):

        count_variants_ena = self.queries.count_variants(cache=False, source=DataSource.ENA.name)
        count_variants_gisaid = self.queries.count_variants(cache=False, source=DataSource.GISAID.name)
        count_samples_ena = self.queries.count_samples(source=DataSource.ENA.name, cache=False)
        count_samples_gisaid = self.queries.count_samples(source=DataSource.GISAID.name, cache=False)
        count_variant_observations_ena = self.queries.count_variant_observations(
            source=DataSource.ENA.name, cache=False)
        count_variant_observations_gisaid = self.queries.count_variant_observations(
            source=DataSource.GISAID.name, cache=False)
        count_subclonal_variant_observations = self.queries.count_subclonal_variant_observations(cache=False)
        count_subclonal_variant_unique = self.queries.count_unique_subclonal_variant(cache=False)
        count_subclonal_variant_unique_only_subclonal = self.queries.count_unique_only_subclonal_variant(cache=False)
        count_countries_ena = self.queries.count_countries(source=DataSource.ENA.name, cache=False)
        count_countries_gisaid = self.queries.count_countries(source=DataSource.GISAID.name, cache=False)

        # delete all rows before starting
        self.session.query(PrecomputedTableCounts).delete()
        self.session.commit()

        database_rows = [
            PrecomputedTableCounts(
                table=Variant.__name__, count=count_variants_ena,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(
                table=GisaidVariant.__name__, count=count_variants_gisaid,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
            PrecomputedTableCounts(
                table=VariantObservation.__name__, count=count_variant_observations_ena,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(
                table=GisaidVariantObservation.__name__, count=count_variant_observations_gisaid,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__, count=count_subclonal_variant_observations),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__ + "_unique", count=count_subclonal_variant_unique),
            PrecomputedTableCounts(
                table=SubclonalVariantObservation.__name__ + "_unique_only_subclonal",
                count=count_subclonal_variant_unique_only_subclonal),
            PrecomputedTableCounts(
                table=SampleEna.__name__, count=count_samples_ena,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(
                table=SampleGisaid.__name__, count=count_samples_gisaid,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
            PrecomputedTableCounts(
                table=PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY, count=count_countries_ena,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.ENA.name),
            PrecomputedTableCounts(
                table=PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY, count=count_countries_gisaid,
                factor=PrecomputedTableCounts.FACTOR_SOURCE, value=DataSource.GISAID.name),
        ]

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedTableCounts.__tablename__))

    def load_variant_abundance_histogram(self):

        database_rows = self._read_variant_abundance_histogram()

        # delete all rows before starting
        self.session.query(PrecomputedVariantAbundanceHistogram).delete()
        self.session.commit()

        if len(database_rows) > 0:
            self.session.add_all(database_rows)
            self.session.commit()
        logger.info("Added {} entries to {}".format(
            len(database_rows), PrecomputedVariantAbundanceHistogram.__tablename__))

    def _read_variant_abundance_histogram(self):
        histograms = []
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
        return database_rows
