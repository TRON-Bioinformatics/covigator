from datetime import date, datetime
from typing import List, Union
import pandas as pd
from scipy.spatial.distance import squareform
from logzero import logger
from sklearn import manifold
from sklearn.cluster import DBSCAN
from sqlalchemy import and_, desc, asc, func, or_, String, text, DateTime, cast
from sqlalchemy.engine.default import DefaultDialect
from sqlalchemy.orm import Session, aliased
from sqlalchemy.sql.sqltypes import NullType, Numeric, Float

from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna, JobEna, JobStatus, VariantObservation, \
    Gene, Variant, VariantCooccurrence, Conservation, JobGisaid

SYNONYMOUS_VARIANT = "synonymous_variant"


class Queries:

    FINAL_JOB_STATE = JobStatus.FINISHED

    def __init__(self, session: Session):
        self.session = session

    def find_job_by_accession_and_status(
            self, run_accession: str, status: JobStatus, data_source: DataSource) -> Union[JobEna, JobGisaid]:
        if data_source == DataSource.ENA:
            return self.session.query(JobEna) \
                .filter(and_(JobEna.run_accession == run_accession, JobEna.status == status)) \
                .first()
        elif data_source == DataSource.GISAID:
            return self.session.query(JobGisaid)\
                .filter(and_(JobGisaid.run_accession == run_accession, JobGisaid.status == status)) \
                .first()
        else:
            raise ValueError("Bad data source {}".format(data_source))

    def find_job_by_accession(self, run_accession: str, data_source: DataSource) -> Union[JobEna, JobGisaid]:
        if data_source == DataSource.ENA:
            return self.session.query(JobEna).filter(JobEna.run_accession == run_accession).first()
        elif data_source == DataSource.GISAID:
            return self.session.query(JobGisaid).filter(JobGisaid.run_accession == run_accession).first()
        else:
            raise ValueError("Bad data source {}".format(data_source))

    def find_first_pending_job(self, data_source: DataSource) -> Union[JobEna, JobGisaid]:
        if data_source == DataSource.ENA:
            return self.session.query(JobEna) \
                .filter(JobEna.status == JobStatus.PENDING) \
                .order_by(JobEna.created_at.desc()) \
                .first()
        elif data_source == DataSource.GISAID:
            return self.session.query(JobGisaid) \
                .filter(JobGisaid.status == JobStatus.PENDING) \
                .order_by(JobGisaid.created_at.desc()) \
                .first()
        else:
            raise ValueError("Bad data source {}".format(data_source))

    def count_jobs_in_queue(self, data_source):
        if data_source == DataSource.ENA:
            count = self.session.query(JobEna).filter(JobEna.status.in_((
                JobStatus.QUEUED, JobStatus.DOWNLOADED, JobStatus.PROCESSED, JobStatus.LOADED))).count()
        elif data_source == DataSource.GISAID:
            count = self.session.query(JobGisaid).filter(JobGisaid.status.in_((
                JobStatus.QUEUED, JobStatus.DOWNLOADED, JobStatus.PROCESSED, JobStatus.LOADED))).count()
        else:
            raise ValueError("Bad data source {}".format(data_source))
        return count

    def find_sample_ena_by_accession(self, run_accession: str) -> SampleEna:
        return self.session.query(SampleEna).filter(SampleEna.run_accession == run_accession).first()

    def get_accumulated_samples_by_country(self) -> pd.DataFrame:
        """
        Returns a DataFrame with columns: data, country, cumsum, count
        """
        samples = pd.read_sql(self.session.query(SampleEna).join(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE).statement,
                              self.session.bind)

        filled_table = None
        if samples.shape[0] > 0:
            # merge countries with less than 10 samples into OTHER
            country_value_counts = samples.country.value_counts()
            other_countries = list(country_value_counts[country_value_counts < 10].index)
            null_country_values = ["none", "not available"]
            samples["country_merged"] = samples.country.transform(
                lambda c: "Other" if c in other_countries or c is None or c.lower() in null_country_values else c)

            # counts samples by country and data
            sample_counts = samples[["first_created", "run_accession", "country_merged"]] \
                .groupby(["first_created", "country_merged"]).count()
            sample_counts.reset_index(inplace=True)
            sample_counts.rename(columns={"run_accession": "count"}, inplace=True)

            # accumulates count ordered by date
            sample_counts['cumsum'] = sample_counts.groupby(['country_merged'])['count'].cumsum()

            # creates empty table with all pairwise combinations of date and country
            dates = sample_counts.first_created.unique()
            countries = sample_counts.sort_values("cumsum", ascending=False).country_merged.unique()
            empty_table = pd.DataFrame(
                index=pd.MultiIndex.from_product([dates, countries], names=["first_created", "country_merged"]))
            empty_table["count"] = 0

            # adds values into empty table
            filled_table = empty_table + sample_counts.set_index(["first_created", "country_merged"])
            filled_table.fillna(0, inplace=True)
            filled_table.reset_index(inplace=True)
            filled_table['cumsum'] = filled_table.groupby(['country_merged'])['count'].cumsum()
            filled_table.rename(columns={"first_created": "date", "country_merged": "country"}, inplace=True)

        return filled_table

    def get_sample_months(self, pattern) -> List[datetime]:
        dates = [
            d[0] for d in
            self.session.query(SampleEna.first_created).join(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE).all()]
        return sorted(set([d.strftime(pattern) for d in dates]))

    def get_gene(self, gene_name: str):
        return self.session.query(Gene).filter(Gene.name == gene_name).first()

    def get_genes(self):
        return [g[0] for g in self.session.query(Gene.name).order_by(Gene.start).all()]

    def get_genes_metadata(self):
        return self.session.query(Gene).order_by(Gene.start).all()

    def get_non_synonymous_variants_by_region(self, start, end) -> pd.DataFrame:
        subquery = self.session.query(VariantObservation.position, Variant.annotation, Variant.hgvs_p,
                               func.count(VariantObservation.position).label("count_occurrences"))\
            .join(Variant)\
            .filter(and_(Variant.position >= start, Variant.position <= end,
                         Variant.annotation != SYNONYMOUS_VARIANT))\
            .group_by(VariantObservation.position, Variant.annotation, Variant.hgvs_p).subquery()
        return pd.read_sql(
            self.session.query(subquery).filter(subquery.c.count_occurrences > 1).statement, self.session.bind)

    def get_variants_by_sample(self, sample_id) -> List[VariantObservation]:
        return self.session.query(VariantObservation) \
            .filter(VariantObservation.sample == sample_id) \
            .order_by(VariantObservation.position, VariantObservation.reference, VariantObservation.alternate) \
            .all()

    def get_variant_cooccurrence(self, variant_one: Variant, variant_two: Variant) -> VariantCooccurrence:
        return self.session.query(VariantCooccurrence) \
            .filter(and_(VariantCooccurrence.variant_id_one == variant_one.variant_id,
                         VariantCooccurrence.variant_id_two == variant_two.variant_id)) \
            .first()

    def count_ena_samples(self) -> int:
        return self.session.query(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE).count()

    def count_countries(self):
        return self.session.query(SampleEna).join(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE)\
            .distinct(SampleEna.country).count()

    def count_variants(self):
        return self.session.query(Variant).count()

    def count_variant_observations(self):
        return self.session.query(VariantObservation).count()

    def get_date_of_first_ena_sample(self) -> date:
        """
        Returns the date of the earliest ENA sample loaded in the database
        """
        result = self.session.query(SampleEna.first_created).join(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE) \
            .order_by(asc(SampleEna.first_created)).first()
        return result[0] if result is not None else result

    def get_date_of_most_recent_ena_sample(self) -> date:
        """
        Returns the date of the latest ENA sample loaded in the database
        """
        result = self.session.query(SampleEna.first_created).join(JobEna).filter(JobEna.status == self.FINAL_JOB_STATE) \
            .order_by(desc(SampleEna.first_created)).first()
        return result[0] if result is not None else result

    def get_date_of_last_check(self, data_source: DataSource) -> date:
        """
        Returns the date of the latest non failed accessor check that also has a subsequent non failed processor run.
        Until the processor has ran the data fetched from the accessor is not available.
        """
        result1 = self.session.query(Log.start).filter(
            and_(Log.source == data_source, Log.module == CovigatorModule.PROCESSOR, Log.has_error == False)).order_by(
            desc(Log.start)).first()
        most_recent_processor_run = result1[0] if result1 is not None else result1
        result2 = None
        if most_recent_processor_run:
            result2 = self.session.query(Log.start).filter(
                and_(Log.source == data_source, Log.module == CovigatorModule.ACCESSOR, Log.has_error == False,
                     Log.start < most_recent_processor_run)).order_by(desc(Log.start)).first()
        return result2[0] if result2 is not None else result2

    def get_date_of_last_update(self, data_source: DataSource) -> date:
        """
        Returns the date of the latest non failed accessor check **with some new data** that also has a subsequent non
        failed processor run.
        Until the processor has ran the data fetched from the accessor is not available.
        """
        result1 = self.session.query(Log.start).filter(
            and_(Log.source == data_source, Log.module == CovigatorModule.PROCESSOR, Log.has_error == False)).order_by(
            desc(Log.start)).first()
        most_recent_processor_run = result1[0] if result1 is not None else result1
        result2 = None
        if most_recent_processor_run:
            result2 = self.session.query(Log.start).filter(
                and_(Log.source == data_source, Log.module == CovigatorModule.ACCESSOR, Log.has_error == False,
                     Log.processed > 0, Log.start < most_recent_processor_run)) \
                .order_by(desc(Log.start)).first()
        return result2[0] if result2 is not None else result2

    def get_top_occurring_variants(self, top=10, gene_name=None, metric="count") -> pd.DataFrame:
        """
        Returns the top occurring variants + the segregated counts of occurrences per month
        with columns: chromosome, position, reference, alternate, total, month, count
        """
        # query for top occurring variants
        query = self.session.query(
            VariantObservation.variant_id, Variant.hgvs_p, Variant.gene_name, Variant.annotation,
            func.count().label('total')).join(Variant)
        if gene_name is not None:
            query = query.filter(and_(Variant.gene_name == gene_name, Variant.annotation != SYNONYMOUS_VARIANT))
        else:
            query = query.filter(Variant.annotation != SYNONYMOUS_VARIANT)
        query = query.group_by(
            VariantObservation.variant_id, Variant.hgvs_p, Variant.gene_name, Variant.annotation)
        query = query.order_by(desc('total')).limit(top)
        top_occurring_variants = pd.read_sql(query.statement, self.session.bind)

        if top_occurring_variants is not None and top_occurring_variants.shape[0] > 0:
            variant_counts_by_month = []

            # NOTE: one query per variant for the counts per month, will this be efficient?
            for _, variant in top_occurring_variants.iterrows():
                variant_counts_by_month.append(self.get_variant_counts_by_month(variant.variant_id))
            top_occurring_variants_by_month = pd.concat(variant_counts_by_month)

            # get total count of samples per month to calculate the frequency by month
            sample_counts_by_month = self.get_sample_counts_by_month()
            top_occurring_variants_by_month = pd.merge(
                left=top_occurring_variants_by_month, right=sample_counts_by_month, how="left", on="month")
            top_occurring_variants_by_month["frequency_by_month"] = \
                (top_occurring_variants_by_month["count"] / top_occurring_variants_by_month["sample_count"]).\
                transform(lambda x: round(x, 3))

            # join both tables with total counts and counts per month
            top_occurring_variants = pd.merge(
                left=top_occurring_variants, right=top_occurring_variants_by_month, on="variant_id", how="left")

            # format the month column appropriately
            top_occurring_variants.month = top_occurring_variants.month.transform(
                lambda d: "{}-{:02d}".format(d.year, int(d.month)))

            # formats the DNA mutation
            top_occurring_variants.rename(columns={'variant_id': 'dna_mutation'}, inplace=True)

            # replace the total count by the frequency
            count_samples = self.count_ena_samples()
            top_occurring_variants['frequency'] = top_occurring_variants.total.transform(
                lambda t: round(float(t) / count_samples, 3))

            # pivots the table over months
            top_occurring_variants = pd.pivot_table(
                top_occurring_variants, index=['gene_name', 'dna_mutation', 'hgvs_p', 'annotation', "frequency", "total"],
                columns=["month"], values=[metric], fill_value=0).droplevel(0, axis=1).reset_index()

        return top_occurring_variants

    def get_variant_counts_by_month(self, variant_id) -> pd.DataFrame:
        query = self.session.query(
            VariantObservation.variant_id,
            func.date_trunc('month', SampleEna.first_created).label("month"),
            func.count().label("count"))\
            .join(SampleEna, VariantObservation.sample == SampleEna.run_accession) \
            .join(JobEna, VariantObservation.sample == JobEna.run_accession) \
            .filter(and_(VariantObservation.variant_id == variant_id, JobEna.status == JobStatus.FINISHED))\
            .group_by(VariantObservation.variant_id, func.date_trunc('month', SampleEna.first_created))
        return pd.read_sql(query.statement, self.session.bind)

    def get_sample_counts_by_month(self) -> pd.DataFrame:
        query = self.session.query(
            func.date_trunc('month', SampleEna.first_created).label("month"),
            func.count().label("sample_count"))\
            .join(JobEna) \
            .filter(JobEna.status == JobStatus.FINISHED) \
            .group_by(func.date_trunc('month', SampleEna.first_created))
        return pd.read_sql(query.statement, self.session.bind)

    def get_variants_cooccurrence_by_gene(self, gene_name, min_cooccurrence=5, test=False) -> pd.DataFrame:
        """
        Returns the full cooccurrence matrix of all non synonymous variants in a gene with at least
        min_occurrences occurrences.
        """
        # query for total samples required to calculate frequencies
        count_samples = self.count_ena_samples()

        # query for cooccurrence matrix
        variant_one = aliased(Variant)
        variant_two = aliased(Variant)
        if not test:
            query = self.session.query(VariantCooccurrence,
                                       variant_one.position,
                                       variant_one.reference,
                                       variant_one.alternate,
                                       variant_one.hgvs_p,
                                       (variant_one.hgvs_p + " - " + variant_two.hgvs_p).label("hgvs_tooltip"))
        else:
            # this is needed for testing environment with SQLite
            query = self.session.query(VariantCooccurrence,
                                       variant_one.position,
                                       variant_one.reference,
                                       variant_one.alternate,
                                       variant_one.hgvs_p,
                                       (variant_one.hgvs_p + " - " + variant_two.hgvs_p).label("hgvs_tooltip"))

        query = query.filter(VariantCooccurrence.count >= min_cooccurrence)\
            .join(variant_one, and_(VariantCooccurrence.variant_id_one == variant_one.variant_id)) \
            .join(variant_two, and_(VariantCooccurrence.variant_id_two == variant_two.variant_id))

        if gene_name is not None:
            query = query.filter(and_(variant_one.gene_name == gene_name,
                         variant_one.annotation != SYNONYMOUS_VARIANT,
                         variant_two.gene_name == gene_name,
                         variant_two.annotation != SYNONYMOUS_VARIANT))
        else:
            query = query.filter(and_(variant_one.annotation != SYNONYMOUS_VARIANT,
                                      variant_two.annotation != SYNONYMOUS_VARIANT))

        self._print_query(query=query)
        data = pd.read_sql(query.statement, self.session.bind)

        full_matrix = None
        if data.shape[0] > 0:
            # these are views of the original data
            annotations = data.loc[data.variant_id_one == data.variant_id_two,
                                   ["variant_id_one", "position", "reference", "alternate", "hgvs_p"]]
            tooltip = data.loc[:, ["variant_id_one", "variant_id_two", "hgvs_tooltip"]]
            diagonal = data.loc[data.variant_id_one == data.variant_id_two, ["variant_id_one", "variant_id_two", "count"]]
            sparse_matrix = data.loc[:, ["variant_id_one", "variant_id_two", "count"]]
            sparse_matrix["frequency"] = sparse_matrix["count"] / count_samples
            sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_one", how="left", suffixes=("", "_one"))
            sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_two", how="left", suffixes=("", "_two"))
            sparse_matrix["count_union"] = sparse_matrix["count_one"] + sparse_matrix["count_two"] - sparse_matrix["count"]
            sparse_matrix["jaccard"] = sparse_matrix["count"] / sparse_matrix["count_union"]
            del sparse_matrix["count_union"]
            del sparse_matrix["count_one"]
            del sparse_matrix["count_two"]

            # from the sparse matrix builds in memory the complete matrix
            all_variants = data.variant_id_one.unique()
            empty_full_matrix = pd.DataFrame(index=pd.MultiIndex.from_product(
                [all_variants, all_variants], names=["variant_id_one", "variant_id_two"])).reset_index()
            full_matrix = pd.merge(
                left=empty_full_matrix, right=sparse_matrix, on=["variant_id_one", "variant_id_two"], how='left')

            # add annotation on variant one
            full_matrix = pd.merge(left=full_matrix, right=annotations, on="variant_id_one", how='left')\
                .rename(columns={"position": "position_one",
                                 "reference": "reference_one",
                                 "alternate": "alternate_one",
                                 "hgvs_p": "hgvs_p_one"})

            # add annotations on variant two
            full_matrix = pd.merge(left=full_matrix, right=annotations,
                                   left_on="variant_id_two", right_on="variant_id_one", how='left') \
                .rename(columns={"variant_id_one_x": "variant_id_one",
                                 "position": "position_two",
                                 "reference": "reference_two",
                                 "alternate": "alternate_two",
                                 "hgvs_p": "hgvs_p_two"})

            # add tooltip
            full_matrix = pd.merge(left=full_matrix, right=tooltip, on=["variant_id_one", "variant_id_two"], how='left')
            # correct tooltip in diagonal
            full_matrix.loc[full_matrix.variant_id_one == full_matrix.variant_id_two, "hgvs_tooltip"] = full_matrix.hgvs_p_one

            # NOTE: transpose matrix manually as plotly transpose does not work with labels
            # the database return the upper diagonal, the lower is best for plots
            full_matrix.sort_values(["position_two", "reference_two", "alternate_two",
                                     "position_one", "reference_one", "alternate_one"], inplace=True)

            full_matrix = full_matrix.loc[:, ["variant_id_one", "variant_id_two", "count", "frequency", "jaccard", "hgvs_tooltip"]]

        return full_matrix

    def get_mds(self, gene_name) -> pd.DataFrame:
        variant_one = aliased(Variant)
        variant_two = aliased(Variant)
        query = self.session.query(VariantCooccurrence)\
            .join(variant_one, and_(VariantCooccurrence.variant_id_one == variant_one.variant_id)) \
            .join(variant_two, and_(VariantCooccurrence.variant_id_two == variant_two.variant_id))

        if gene_name is not None:
            query = query.filter(and_(variant_one.gene_name == gene_name,
                                      variant_one.annotation != SYNONYMOUS_VARIANT,
                                      variant_two.gene_name == gene_name,
                                      variant_two.annotation != SYNONYMOUS_VARIANT))
        else:
            query = query.filter(and_(variant_one.annotation != SYNONYMOUS_VARIANT,
                                      variant_two.annotation != SYNONYMOUS_VARIANT))

        sparse_matrix = pd.read_sql(query.statement, self.session.bind)
        diagonal = sparse_matrix.loc[sparse_matrix.variant_id_one == sparse_matrix.variant_id_two,
                                     ["variant_id_one", "variant_id_two", "count"]]
        sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_one", how="left",
                                 suffixes=("", "_one"))
        sparse_matrix = pd.merge(left=sparse_matrix, right=diagonal, on="variant_id_two", how="left",
                                 suffixes=("", "_two"))
        sparse_matrix["count_union"] = sparse_matrix["count_one"] + sparse_matrix["count_two"] - sparse_matrix["count"]
        sparse_matrix["jaccard_dissimilarity"] = 1 - sparse_matrix["count"] / sparse_matrix["count_union"]

        all_variants = sparse_matrix.variant_id_one.unique()
        empty_full_matrix = pd.DataFrame(index=pd.MultiIndex.from_product(
            [all_variants, all_variants], names=["variant_id_one", "variant_id_two"])).reset_index()
        full_matrix = pd.merge(
            # gets only the inferior matrix without the diagnonal
            left=empty_full_matrix.loc[empty_full_matrix.variant_id_one < empty_full_matrix.variant_id_two, :],
            right=sparse_matrix.loc[:, ["variant_id_one", "variant_id_two", "jaccard_dissimilarity"]],
            on=["variant_id_one", "variant_id_two"], how='left')
        full_matrix.fillna(1.0, inplace=True)
        # TODO: sort this as described here https://stackoverflow.com/questions/13079563/how-does-condensed-distance-matrix-work-pdist
        full_matrix.sort_values(by=["variant_id_one", "variant_id_two"], inplace=True)

        # TODO_ replace this call by http://scikit-bio.org/docs/0.4.2/generated/generated/skbio.stats.distance.DissimilarityMatrix.html
        # to maintain variant labels
        distance_matrix = squareform(full_matrix.jaccard_dissimilarity)
        mds_model = manifold.MDS(n_components=3, random_state=123, dissimilarity='precomputed')
        mds_coords = mds_model.fit_transform(distance_matrix)

        # performs clustering
        clusters = DBSCAN().fit_predict(distance_matrix)

        # builds data into a dataframe
        data = pd.DataFrame(mds_coords, columns=["PC1", "PC2", "PC3"])
        data["cluster"] = clusters

        return data

    def get_variant_abundance_histogram(self, bin_size=50) -> pd.DataFrame:

        # queries for the maximum position
        maximum_position = self.session.query(func.max(Variant.position)).first()[0]
        histogram = None
        if maximum_position is not None:
            # builds all possible bins
            all_bins = pd.DataFrame(data=[i*bin_size for i in range(int(maximum_position/bin_size) + 1)], columns=["position_bin"])

            # counts variants over those bins
            sql_query = """
                    SELECT cast("position"/{bin_size} as int)*{bin_size} AS position_bin,
                           COUNT(*) as count_unique_variants
                    FROM {table_name}
                    GROUP BY position_bin
                    ORDER BY position_bin;
                    """.format(bin_size=bin_size, table_name=Variant.__tablename__)
            binned_counts_variants = pd.read_sql_query(sql_query, self.session.bind)

            # counts variant observations over those bins
            sql_query = """
                    SELECT cast("position"/{bin_size} as int)*{bin_size} AS position_bin,
                           COUNT(*) as count_variant_observations
                    FROM {table_name}
                    GROUP BY position_bin
                    ORDER BY position_bin;
                    """.format(bin_size=bin_size, table_name=VariantObservation.__tablename__)
            binned_counts_variant_observations = pd.read_sql_query(sql_query, self.session.bind)

            histogram = all_bins.set_index("position_bin").join(binned_counts_variants.set_index("position_bin"))
            histogram = histogram.join(binned_counts_variant_observations.set_index("position_bin"))
            histogram.fillna(0, inplace=True)
            histogram.reset_index(inplace=True)

        return histogram

    def get_conservation_table(self, bin_size=50, start=None, end=None) -> pd.DataFrame:
        # counts variants over those bins
        sql_query = """
                SELECT cast("start"/{bin_size} as int)*{bin_size} AS position_bin,
                       AVG("conservation") as conservation,
                       AVG("conservation_sarbecovirus") as conservation_sarbecovirus,
                       AVG("conservation_vertebrates") as conservation_vertebrates
                FROM {table_name}
                {where}
                GROUP BY position_bin
                ORDER BY position_bin;
                """.format(bin_size=bin_size, table_name= Conservation.__tablename__,
                           where="WHERE start >= {start} and start <= {end}".format(start=start, end=end)
                           if start is not None and end is not None else "")
        return pd.read_sql_query(sql_query, self.session.bind)

    def _print_query(self, query):
        class StringLiteral(String):
            """Teach SA how to literalize various things."""

            def literal_processor(self, dialect):
                super_processor = super(StringLiteral, self).literal_processor(dialect)

                def process(value):
                    if isinstance(value, int):
                        return str(value)
                    if not isinstance(value, str):
                        value = str(value)
                    result = super_processor(value)
                    if isinstance(result, bytes):
                        result = result.decode(dialect.encoding)
                    return result

                return process

        class LiteralDialect(DefaultDialect):
            colspecs = {
                # prevent various encoding explosions
                String: StringLiteral,
                # teach SA about how to literalize a datetime
                DateTime: StringLiteral,
                # don't format py2 long integers to NULL
                NullType: StringLiteral,
            }
        logger.info(query.statement.compile(
            dialect=LiteralDialect(),
            compile_kwargs={'literal_binds': True}).string)
