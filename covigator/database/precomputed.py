import pandas as pd
from sqlalchemy.orm import Session
from covigator.configuration import Configuration
from covigator.database.database import get_database
from covigator.database.model import PrecomputedVariantsPerSample


class Precomputer:

    def __init__(self, session: Session):
        self.session = session

    def load_counts_variants_per_sample(self):

        # reads the data from the database
        sql_query = """
        select count(*), counts.number_mutations, counts.source, counts.variant_type, counts.gene_name from (
            select count(*) as number_mutations, source, variant_type, gene_name, sample from variant_observation_v13
            group by source, variant_type, gene_name, sample)
            as counts group by counts.number_mutations, counts.source, counts.variant_type, counts.gene_name;
            """
        data = pd.read_sql_query(sql_query, self.session.bind)
        sql_query = """
        select count(*), counts.number_mutations, counts.source, counts.variant_type from (
            select count(*) as number_mutations, source, variant_type, sample from variant_observation_v13
            group by source, variant_type, sample)
            as counts group by counts.number_mutations, counts.source, counts.variant_type;
            """
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


if __name__ == '__main__':
    database = get_database(config=Configuration(), initialize=True, verbose=True)
    Precomputer(session=database.get_database_session()).load_counts_variants_per_sample()
