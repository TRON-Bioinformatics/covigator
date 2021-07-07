import pandas as pd
from sqlalchemy import func, desc, and_
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.database.database import get_database, Database
from covigator.database.model import PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, \
    PRECOMPUTED_VARIANTS_PER_SAMPLE_TABLE_NAME, VARIANT_OBSERVATION_TABLE_NAME, PrecomputedIndelLength, \
    PrecomputedAnnotation, VariantObservation, DataSource, PrecomputedOccurrence
from logzero import logger

from covigator.database.queries import SYNONYMOUS_VARIANT, Queries

NUMBER_TOP_OCCURRENCES = 1000


class Precomputer:

    def __init__(self, session: Session):
        self.session = session

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
        queries = Queries(self.session)

        # gets the top occurrent variants for each source and overall
        top_occurring_variants_ena = queries.get_top_occurring_variants(
            top=NUMBER_TOP_OCCURRENCES, source=DataSource.ENA)
        top_occurring_variants_gisaid = queries.get_top_occurring_variants(
            top=NUMBER_TOP_OCCURRENCES, source=DataSource.GISAID)
        top_occurring_variants = queries.get_top_occurring_variants(top=NUMBER_TOP_OCCURRENCES, source=None)

        # delete all rows before starting
        self.session.query(PrecomputedOccurrence).delete()
        self.session.commit()

        # stores the precomputed data
        database_rows = []
        for index, row in top_occurring_variants_ena.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedOccurrence(
                total=row["total"],
                frequency=row["frequency"],
                variant_id=row["variant_id"],
                hgvs_p=row["hgvs_p"],
                gene_name=row["gene_name"],
                annotation=row["annotation_highest_impact"],
                source=DataSource.ENA,
                month=row["month"],
                count=row["count"],
                frequency_by_month=row["frequency_by_month"],
            ))

        for index, row in top_occurring_variants_gisaid.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedOccurrence(
                total=row["total"],
                frequency=row["frequency"],
                variant_id=row["variant_id"],
                hgvs_p=row["hgvs_p"],
                gene_name=row["gene_name"],
                annotation=row["annotation_highest_impact"],
                source=DataSource.GISAID,
                month=row["month"],
                count=row["count"],
                frequency_by_month=row["frequency_by_month"],
            ))

        for index, row in top_occurring_variants.iterrows():
            # add entries per gene
            database_rows.append(PrecomputedOccurrence(
                total=row["total"],
                frequency=row["frequency"],
                variant_id=row["variant_id"],
                hgvs_p=row["hgvs_p"],
                gene_name=row["gene_name"],
                annotation=row["annotation_highest_impact"],
                month=row["month"],
                count=row["count"],
                frequency_by_month=row["frequency_by_month"],
            ))

        self.session.add_all(database_rows)
        self.session.commit()
        logger.info("Added {} entries to {}".format(len(database_rows), PrecomputedOccurrence.__tablename__))


if __name__ == '__main__':
    database = Database(initialize=True, config=Configuration())
    Precomputer(session=database.get_database_session()).load_top_occurrences()
