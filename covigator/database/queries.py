from datetime import date, datetime
from typing import List, Union
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from logzero import logger
from skbio.stats.distance import DissimilarityMatrix
from sklearn.cluster import OPTICS
from sklearn.manifold import MDS
from sqlalchemy import and_, desc, asc, func, String, DateTime, cast
from sqlalchemy.engine.default import DefaultDialect
from sqlalchemy.orm import Session, aliased
from sqlalchemy.sql.sqltypes import NullType

from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna, JobEna, JobStatus, VariantObservation, \
    Gene, Variant, VariantCooccurrence, Conservation, JobGisaid, SampleGisaid, SubclonalVariantObservation, \
    PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, PrecomputedIndelLength, VariantType, \
    PrecomputedAnnotation
from covigator.exceptions import CovigatorQueryException

SYNONYMOUS_VARIANT = "synonymous_variant"


class Queries:

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

    def find_sample_by_accession(self, run_accession: str, source: DataSource) -> Union[SampleEna, SampleGisaid]:
        if source == DataSource.ENA:
            sample = self.session.query(SampleEna).filter(SampleEna.run_accession == run_accession).first()
        elif source == DataSource.GISAID:
            sample = self.session.query(SampleGisaid).filter(SampleGisaid.run_accession == run_accession).first()
        else:
            raise CovigatorQueryException("Bad query trying to fetch a sample")
        return sample

    def get_countries(self) -> List[str]:
        countries_ena = [c for c, in self.session.query(SampleEna.country).filter(
            SampleEna.finished).distinct().all()]
        countries_gisaid = [c for c, in self.session.query(SampleGisaid.country).filter(
            SampleGisaid.finished).distinct().all()]
        return sorted(list(set(countries_ena + countries_gisaid)))

    def get_variants_per_sample(self, data_source: str, genes: List[str]):
        """
        Returns a DataFrame with columns: number_mutations, count, type
        where type: SNV, insertion or deletion
        """
        query = self.session.query(PrecomputedVariantsPerSample)
        if data_source is not None:
            query = query.filter(PrecomputedVariantsPerSample.source == data_source)
        if genes is not None and genes:
            query = query.filter(PrecomputedVariantsPerSample.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedVariantsPerSample.gene_name == None)

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            data.variant_type = data.variant_type.transform(lambda x: x.name if x else x)
        return data[["number_mutations", "variant_type", "count"]] \
            .groupby(["number_mutations", "variant_type"]).sum().reset_index()

    def get_indel_lengths(self, data_source, genes):
        query = self.session.query(PrecomputedIndelLength)
        if data_source is not None:
            query = query.filter(PrecomputedIndelLength.source == data_source)
        if genes is not None and genes:
            query = query.filter(PrecomputedIndelLength.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedIndelLength.gene_name == None)

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            data["tmp_variant_type"] = data["length"].transform(
                lambda x: VariantType.INSERTION.name if x > 0 else VariantType.DELETION.name)
            data["inframe"] = data["length"].transform(lambda x: "INFRAME" if x % 3 == 0 else "FRAMESHIFT")
            data["variant_type"] = data[["tmp_variant_type", "inframe"]].apply(
                lambda x: "{}_{}".format(x[0], x[1]), axis=1)
            return data[["length", "variant_type", "count"]] \
                .groupby(["length", "variant_type", ]).sum().reset_index() \
                .sort_values("length", ascending=True)
        return data

    def get_annotations(self, data_source, genes):
        query = self.session.query(PrecomputedAnnotation)
        if data_source is not None:
            query = query.filter(PrecomputedAnnotation.source == data_source)
        if genes is not None and genes:
            query = query.filter(PrecomputedAnnotation.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedAnnotation.gene_name == None)

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            return data[["annotation", "count"]] \
                .groupby(["annotation", ]).sum().reset_index() \
                .sort_values("count", ascending=False)
        return data

    def get_substitutions(self, data_source, genes, variant_types):
        query = self.session.query(PrecomputedSubstitutionsCounts)
        if variant_types is not None and variant_types:
            query = query.filter(PrecomputedSubstitutionsCounts.variant_type.in_(variant_types))
        if data_source is not None:
            query = query.filter(PrecomputedSubstitutionsCounts.source == data_source)
        if genes is not None and genes:
            query = query.filter(PrecomputedSubstitutionsCounts.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedSubstitutionsCounts.gene_name == None)

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            data.variant_type = data.variant_type.transform(lambda x: x.name if x else x)
            data["substitution"] = data[["reference", "alternate"]].apply(lambda x: "{}>{}".format(x[0], x[1]), axis=1)
            data = data[["substitution", "variant_type", "count"]] \
                .groupby(["substitution", "variant_type"]).sum().reset_index()\
                .sort_values("count", ascending=False)
            data["rate"] = (data["count"] / data["count"].sum()).transform(lambda x: "{} %".format(round(x * 100, 1)))
            data = data.head(20)
        return data

    def get_accumulated_samples_by_country(
            self, data_source: DataSource, countries: List[str], min_samples=100) -> pd.DataFrame:
        """
        Returns a DataFrame with columns: data, country, cumsum, count
        """
        samples_ena = None

        if data_source is None or data_source == DataSource.ENA.name:
            query = self.session.query(
                func.count().label("count"), SampleEna.first_created.label("date"), SampleEna.country) \
                .filter(SampleEna.finished) \
                .group_by(SampleEna.first_created, SampleEna.country)
            if countries:
                query = query.filter(SampleEna.country.in_(countries))
            samples_ena = pd.read_sql(query.statement, self.session.bind).astype(
                {'date': 'datetime64', 'count': 'float64'})

        samples_gisaid = None
        if data_source is None or data_source == DataSource.GISAID.name:
            query = self.session.query(
                func.count().label("count"), SampleGisaid.date, SampleGisaid.country) \
                .filter(SampleGisaid.finished)\
                .group_by(SampleGisaid.date, SampleGisaid.country)
            if countries:
                query = query.filter(SampleGisaid.country.in_(countries))
            samples_gisaid = pd.read_sql(query.statement, self.session.bind).astype(
                {'date': 'datetime64', 'count': 'float64'})

        if samples_gisaid is None:
            samples = samples_ena
        elif samples_ena is None:
            samples = samples_gisaid
        else:
            samples = pd.concat([samples_gisaid, samples_ena]).groupby(['date', 'country']).sum().reset_index()

        filled_table = None
        if samples is not None and samples.shape[0] > 0:
            # accumulates count ordered by date
            samples['cumsum'] = samples.groupby(['country'])['count'].cumsum()

            # creates empty table with all pairwise combinations of date and country
            dates = samples.date[~samples.date.isna()].sort_values().unique()
            countries = samples[(~samples.country.isna()) & (samples["cumsum"] > min_samples)] \
                .sort_values("cumsum", ascending=False).country.unique()

            empty_table = pd.DataFrame(
                index=pd.MultiIndex.from_product([dates, countries], names=["date", "country"]))
            empty_table["count"] = 0

            # adds values into empty table
            filled_table = pd.merge(
                left=empty_table.reset_index(),
                right=samples.loc[:, ["date", "country", "count"]],
                on=["date", "country"],
                how='left',
                suffixes=("_x", "_y")).fillna(0)

            filled_table["count"] = filled_table.count_x + filled_table.count_y
            filled_table['cumsum'] = filled_table.groupby(['country'])['count'].cumsum()

        return filled_table

    def get_sample_months(self, pattern) -> List[datetime]:
        dates_ena = [d.strftime(pattern) for d, in
                     self.session.query(SampleEna.first_created).filter(
                         and_(SampleEna.finished, SampleEna.first_created.isnot(None))).distinct().all()]
        dates_gisaid = [d.strftime(pattern) for d, in
                     self.session.query(SampleGisaid.date).filter(
                         and_(SampleGisaid.finished, SampleGisaid.date.isnot(None))).distinct().all()]
        return sorted(set(dates_ena + dates_gisaid))

    def get_gene(self, gene_name: str):
        return self.session.query(Gene).filter(Gene.name == gene_name).first()

    def get_genes(self):
        return self.session.query(Gene).order_by(Gene.start).all()

    def get_non_synonymous_variants_by_region(self, start, end, source) -> pd.DataFrame:
        query = self.session.query(VariantObservation.position,
                                      VariantObservation.annotation,
                                      VariantObservation.hgvs_p,
                                      func.count().label("count_occurrences"))\
            .filter(and_(VariantObservation.position >= start, VariantObservation.position <= end,
                         VariantObservation.annotation != SYNONYMOUS_VARIANT))
        if source == DataSource.ENA.name:
            query = query.filter(VariantObservation.source == DataSource.ENA)
        elif source == DataSource.GISAID.name:
            query = query.filter(VariantObservation.source == DataSource.GISAID)
        subquery = query.group_by(VariantObservation.position, VariantObservation.annotation, VariantObservation.hgvs_p).subquery()
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

    def count_samples(self, source: str = None) -> int:
        count = 0
        if source is None or source == DataSource.ENA.name:
            count += self.session.query(SampleEna).filter(SampleEna.finished).count()
        if source is None or source == DataSource.GISAID.name:
            count += self.session.query(SampleGisaid).filter(SampleGisaid.finished).count()
        return count

    def count_countries(self):
        return len(self.get_countries())

    def count_variants(self):
        return self.session.query(Variant).count()

    def count_insertions(self):
        return self.session.query(Variant).filter(func.length(Variant.alternate) > 1).count()

    def count_deletions(self):
        return self.session.query(Variant).filter(func.length(Variant.reference) > 1).count()

    def count_variant_observations(self):
        return self.session.query(VariantObservation).count()

    def count_subclonal_variant_observations(self):
        return self.session.query(SubclonalVariantObservation).count()

    def get_date_of_first_sample(self, source: DataSource = DataSource.ENA) -> date:
        """
        Returns the date of the earliest ENA sample loaded in the database
        """
        if source == DataSource.ENA:
            result = self.session.query(SampleEna.first_created).filter(
                and_(SampleEna.finished, SampleEna.first_created.isnot(None))) \
                .order_by(asc(SampleEna.first_created)).first()
        elif source == DataSource.GISAID:
            result = self.session.query(SampleGisaid.date).filter(
                and_(SampleGisaid.finished, SampleGisaid.date.isnot(None))) \
                .order_by(asc(SampleGisaid.date)).first()
        else:
            raise CovigatorQueryException("No valid data source for query of first sample")
        return result[0] if result is not None else result

    def get_date_of_most_recent_sample(self, source: DataSource = DataSource.ENA) -> date:
        """
        Returns the date of the latest ENA sample loaded in the database
        """
        if source == DataSource.ENA:
            result = self.session.query(SampleEna.first_created).filter(
                and_(SampleEna.finished, SampleEna.first_created.isnot(None))) \
                .order_by(desc(SampleEna.first_created)).first()
        elif source == DataSource.GISAID:
            result = self.session.query(SampleGisaid.date).filter(
                and_(SampleGisaid.finished, SampleGisaid.date.isnot(None))) \
                .order_by(desc(SampleGisaid.date)).first()
        else:
            raise CovigatorQueryException("No valid data source for query of most recent sample")
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

    def get_top_occurring_variants(self, top=10, gene_name=None, metric="count", source=None) -> pd.DataFrame:
        """
        Returns the top occurring variants + the segregated counts of occurrences per month
        with columns: chromosome, position, reference, alternate, total, month, count
        """
        # query for top occurring variants
        query = self.session.query(
            VariantObservation.variant_id, VariantObservation.hgvs_p, VariantObservation.gene_name,
            VariantObservation.annotation, func.count().label('total'))
        if gene_name is not None:
            query = query.filter(and_(VariantObservation.gene_name == gene_name,
                                      VariantObservation.annotation != SYNONYMOUS_VARIANT))
        else:
            query = query.filter(VariantObservation.annotation != SYNONYMOUS_VARIANT)
        if source == DataSource.ENA.name:
            query = query.filter(VariantObservation.source == DataSource.ENA)
        elif source == DataSource.GISAID.name:
            query = query.filter(VariantObservation.source == DataSource.GISAID)

        query = query.group_by(
            VariantObservation.variant_id, VariantObservation.hgvs_p, VariantObservation.gene_name,
            VariantObservation.annotation)
        query = query.order_by(desc('total')).limit(top)
        top_occurring_variants = pd.read_sql(query.statement, self.session.bind)

        if top_occurring_variants is not None and top_occurring_variants.shape[0] > 0:
            variant_counts_by_month = []

            # NOTE: one query per variant for the counts per month, will this be efficient?
            for _, variant in top_occurring_variants.iterrows():
                variant_counts_by_month.append(self.get_variant_counts_by_month(variant.variant_id, source=source))
            top_occurring_variants_by_month = pd.concat(variant_counts_by_month)

            # get total count of samples per month to calculate the frequency by month
            sample_counts_by_month = self.get_sample_counts_by_month(source=source)
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
            count_samples = self.count_samples(source=source)
            top_occurring_variants['frequency'] = top_occurring_variants.total.transform(
                lambda t: round(float(t) / count_samples, 3))

            # pivots the table over months
            top_occurring_variants = pd.pivot_table(
                top_occurring_variants, index=['gene_name', 'dna_mutation', 'hgvs_p', 'annotation', "frequency", "total"],
                columns=["month"], values=[metric], fill_value=0).droplevel(0, axis=1).reset_index()

        return top_occurring_variants

    def get_variant_counts_by_month(self, variant_id, source=None) -> pd.DataFrame:
        query = self.session.query(
            VariantObservation.variant_id,
            func.date_trunc('month', VariantObservation.date).label("month"),
            func.count().label("count"))\
            .filter(and_(VariantObservation.variant_id == variant_id, VariantObservation.date.isnot(None)))
        if source == DataSource.ENA.name:
            query = query.filter(VariantObservation.source == DataSource.ENA)
        elif source == DataSource.GISAID.name:
            query = query.filter(VariantObservation.source == DataSource.GISAID)
        query = query.group_by(VariantObservation.variant_id, func.date_trunc('month', VariantObservation.date))
        return pd.read_sql(query.statement, self.session.bind)

    def get_sample_counts_by_month(self, source=None) -> pd.DataFrame:
        counts_ena = None
        if source is None or source == DataSource.ENA.name:
            query = self.session.query(
                func.date_trunc('month', SampleEna.first_created).label("month"),
                func.count().label("sample_count"))\
                .filter(SampleEna.finished) \
                .group_by(func.date_trunc('month', SampleEna.first_created))
            counts_ena = pd.read_sql(query.statement, self.session.bind)
        counts_gisaid = None
        if source is None or source == DataSource.GISAID.name:
            query = self.session.query(
                func.date_trunc('month', SampleGisaid.date).label("month"),
                func.count().label("sample_count"))\
                .filter(SampleGisaid.finished) \
                .group_by(func.date_trunc('month', SampleGisaid.date))
            counts_gisaid = pd.read_sql(query.statement, self.session.bind)
        if counts_gisaid is None:
            # NOTE: setting index and then resetting is necessary to get the date column in the right dtype
            counts = counts_ena
        elif counts_ena is None:
            # NOTE: setting index and then resetting is necessary to get the date column in the right dtype
            counts = counts_gisaid
        else:
            counts = counts_gisaid.set_index(["month"]).add(
                counts_ena.set_index(["month"]), fill_value=0).reset_index()
        return counts

    def get_variants_cooccurrence_by_gene(self, gene_name, min_cooccurrence=5, test=False) -> pd.DataFrame:
        """
        Returns the full cooccurrence matrix of all non synonymous variants in a gene with at least
        min_occurrences occurrences.
        """
        # query for total samples required to calculate frequencies
        # FIXME: once we have the cooccurrence for GISAID samples we need to count all samples here
        count_samples = self.count_samples(source=DataSource.ENA.name)

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

            # calculate Jaccard index
            sparse_matrix["count_union"] = sparse_matrix["count_one"] + sparse_matrix["count_two"] - sparse_matrix["count"]
            sparse_matrix["jaccard"] = sparse_matrix["count"] / sparse_matrix["count_union"]

            # calculate Cohen's kappa
            sparse_matrix["chance_agreement"] = np.exp(-sparse_matrix["count"])
            sparse_matrix["kappa"] = 1 - ((1 - sparse_matrix.jaccard) / (1 - sparse_matrix.chance_agreement))
            sparse_matrix["kappa"] = sparse_matrix["kappa"].transform(lambda k: k if k > 0 else 0)

            del sparse_matrix["count_union"]
            del sparse_matrix["count_one"]
            del sparse_matrix["count_two"]
            del sparse_matrix["chance_agreement"]

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

            full_matrix = full_matrix.loc[:, ["variant_id_one", "variant_id_two", "count", "frequency", "jaccard",
                                              "kappa", "hgvs_tooltip"]]

        return full_matrix

    def get_mds(self, gene_name, min_cooccurrence, min_samples) -> pd.DataFrame:

        variant_one = aliased(Variant)
        variant_two = aliased(Variant)
        query = self.session.query(VariantCooccurrence,
                                   variant_one.hgvs_p.label("hgvs_p_one"),
                                   variant_two.hgvs_p.label("hgvs_p_two")) \
            .filter(VariantCooccurrence.count >= min_cooccurrence) \
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
        sparse_matrix_with_diagonal = pd.merge(
            left=sparse_matrix, right=diagonal, on="variant_id_one", how="left", suffixes=("", "_one"))
        sparse_matrix_with_diagonal = pd.merge(
            left=sparse_matrix_with_diagonal, right=diagonal, on="variant_id_two", how="left", suffixes=("", "_two"))

        # calculate Jaccard index
        sparse_matrix_with_diagonal["count_union"] = sparse_matrix_with_diagonal["count_one"] + \
                                                     sparse_matrix_with_diagonal["count_two"] - \
                                                     sparse_matrix_with_diagonal["count"]
        sparse_matrix_with_diagonal["jaccard_similarity"] = sparse_matrix_with_diagonal["count"] / \
                                                               sparse_matrix_with_diagonal.count_union
        sparse_matrix_with_diagonal["jaccard_dissimilarity"] = 1 - sparse_matrix_with_diagonal.jaccard_similarity

        # calculate Cohen's kappa
        sparse_matrix_with_diagonal["chance_agreement"] = np.exp(-sparse_matrix_with_diagonal["count"])
        sparse_matrix_with_diagonal["kappa"] = 1 - ((1 - sparse_matrix_with_diagonal.jaccard_similarity) / (
                1 - sparse_matrix_with_diagonal.chance_agreement))
        sparse_matrix_with_diagonal["kappa"] = sparse_matrix_with_diagonal["kappa"].transform(lambda k: k if k > 0 else 0)
        sparse_matrix_with_diagonal["kappa_dissimilarity"] = 1 - sparse_matrix_with_diagonal.kappa

        dissimilarity_metric = "kappa_dissimilarity"

        # build upper diagonal matrix
        all_variants = sparse_matrix_with_diagonal.variant_id_one.unique()
        empty_full_matrix = pd.DataFrame(index=pd.MultiIndex.from_product(
            [all_variants, all_variants], names=["variant_id_one", "variant_id_two"])).reset_index()
        upper_diagonal_matrix = pd.merge(
            # gets only the inferior matrix without the diagnonal
            left=empty_full_matrix.loc[empty_full_matrix.variant_id_one < empty_full_matrix.variant_id_two, :],
            right=sparse_matrix_with_diagonal.loc[:, ["variant_id_one", "variant_id_two", dissimilarity_metric]],
            on=["variant_id_one", "variant_id_two"], how='left')
        upper_diagonal_matrix.fillna(1.0, inplace=True)
        upper_diagonal_matrix.sort_values(by=["variant_id_one", "variant_id_two"], inplace=True)

        logger.info("Building square distance matrix...")
        distance_matrix = squareform(upper_diagonal_matrix[dissimilarity_metric])
        # this ensures the order of variants ids is coherent with the non redundant form of the distance matrix
        ids = np.array([list(upper_diagonal_matrix.variant_id_one[0])[0]] + \
              list(upper_diagonal_matrix.variant_id_two[0:len(upper_diagonal_matrix.variant_id_two.unique())]))
        distance_matrix_with_ids = DissimilarityMatrix(data=distance_matrix, ids=ids)

        logger.info("Clustering...")
        clusters = OPTICS(min_samples=min_samples, max_eps=1.4).fit_predict(distance_matrix_with_ids.data)

        logger.info("Dimensionality reduction...")
        dimensionality_reduction_model = MDS(
            n_components=2, random_state=123, dissimilarity='precomputed', n_init=1, max_iter=10)  # this two values make computation faster
        coords = dimensionality_reduction_model.fit_transform(distance_matrix_with_ids.data)

        logger.info("Building clustering dataframe...")
        data = pd.DataFrame(coords, columns=["PC1", "PC2"])
        data["cluster"] = clusters
        data["variant_id"] = distance_matrix_with_ids.ids

        logger.info("Annotate with HGVS.p ...")
        annotations = pd.concat([
            sparse_matrix.loc[:, ["variant_id_one", "hgvs_p_one"]].rename(
                columns={"variant_id_one": "variant_id", "hgvs_p_one": "hgvs_p"}),
            sparse_matrix.loc[:, ["variant_id_two", "hgvs_p_two"]].rename(
                columns={"variant_id_two": "variant_id", "hgvs_p_two": "hgvs_p"})])
        data = pd.merge(left=data, right=annotations, on="variant_id", how="left")
        data["tooltip"] = data[["variant_id", "hgvs_p"]].apply(
            lambda x: "<b>{}</b><br>Variant: {}".format(x[1], x[0]), axis=1)

        logger.info("Annotate with cluster mean Jaccard index...")
        data["cluster_jaccard_mean"] = 1.0
        data["cluster_members"] = 0
        for c in data.cluster.unique():
            variants_in_cluster = data[data.cluster == c].variant_id.unique()
            data.cluster_jaccard_mean = np.where(data.cluster == c, sparse_matrix_with_diagonal[
                (sparse_matrix.variant_id_one.isin(variants_in_cluster)) &
                (sparse_matrix.variant_id_two.isin(variants_in_cluster)) &
                (sparse_matrix.variant_id_one != sparse_matrix.variant_id_two)
            ].jaccard_dissimilarity.mean(), data.cluster_jaccard_mean)
            data.cluster_members = np.where(data.cluster == c, variants_in_cluster.size, data.cluster_members)

        logger.info("Compose tooltip...")
        data["tooltip"] = data[["variant_id", "hgvs_p", "cluster_jaccard_mean", "cluster_members"]].apply(
            lambda x: "<b>{}</b><br>Variant: {}<br>Jaccard mean in cluster:{}<br>Members in cluster: {}".format(
                x[1], x[0], round(1 - x[2], 3), x[3]), axis=1)

        return data

    def get_variant_abundance_histogram(self, bin_size=50, source: str = None) -> pd.DataFrame:

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
                    {source_filter}
                    GROUP BY position_bin
                    ORDER BY position_bin;
                    """.format(bin_size=bin_size, table_name=VariantObservation.__tablename__,
                               source_filter="WHERE source='{}'".format(source) if source is not None else "")
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
