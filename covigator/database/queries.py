from datetime import date, datetime
from typing import List

import pandas as pd
import numpy as np
from sqlalchemy import and_, desc, asc, func
from sqlalchemy.orm import Session, aliased

from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna, JobEna, JobStatus, VariantObservation, \
    Gene, Variant, VariantCooccurrence

SYNONYMOUS_VARIANT = "synonymous_variant"


class Queries:

    FINAL_JOB_STATE = JobStatus.COOCCURRENCE

    def __init__(self, session: Session):
        self.session = session

    def find_job_by_accession_and_status(self, run_accession: str, status: JobStatus) -> JobEna:
        return self.session.query(JobEna)\
            .filter(and_(JobEna.run_accession == run_accession, JobEna.status == status)) \
            .first()

    def find_job_by_accession(self, run_accession: str) -> JobEna:
        return self.session.query(JobEna).filter(JobEna.run_accession == run_accession).first()

    def find_ena_run_by_accession(self, run_accession: str) -> SampleEna:
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
            samples["country_merged"] = samples.country.transform(
                lambda c: "Other" if c in other_countries or c is None or c == "None" else c)

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
        return [g[0] for g in self.session.query(Gene.name).all()]

    def get_non_synonymous_variants_by_gene(self, gene_name) -> pd.DataFrame:
        subquery = self.session.query(VariantObservation.position, Variant.annotation, Variant.hgvs_p,
                               func.count(VariantObservation.position).label("count_occurrences"))\
            .join(Variant)\
            .filter(and_(Variant.gene_name == gene_name, Variant.annotation != SYNONYMOUS_VARIANT))\
            .group_by(VariantObservation.position, Variant.annotation, Variant.hgvs_p).subquery()
        return pd.read_sql(
            self.session.query(subquery).filter(subquery.c.count_occurrences > 1).statement, self.session.bind)

    def get_variants_by_sample(self, sample_id) -> List[VariantObservation]:
        return self.session.query(VariantObservation) \
            .filter(VariantObservation.sample == sample_id) \
            .order_by(VariantObservation.chromosome, VariantObservation.position) \
            .all()

    def get_variant_cooccurrence(self, variant_one: Variant, variant_two: Variant) -> VariantCooccurrence:
        return self.session.query(VariantCooccurrence) \
            .filter(and_(VariantCooccurrence.chromosome_one == variant_one.chromosome,
                         VariantCooccurrence.position_one == variant_one.position,
                         VariantCooccurrence.reference_one == variant_one.reference,
                         VariantCooccurrence.alternate_one == variant_one.alternate,
                         VariantCooccurrence.chromosome_two == variant_two.chromosome,
                         VariantCooccurrence.position_two == variant_two.position,
                         VariantCooccurrence.reference_two == variant_two.reference,
                         VariantCooccurrence.alternate_two == variant_two.alternate)) \
            .first()

    def count_ena_samples_loaded(self) -> int:
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

    def get_top_occurring_variants(self, top=10, gene_name=None) -> pd.DataFrame:
        """
        Returns the top occurring variants + the segregated counts of occurrences per month
        with columns: chromosome, position, reference, alternate, total, month, count
        """
        # query for top occurring variants
        query = self.session.query(
            VariantObservation.position, VariantObservation.reference, VariantObservation.alternate,
            Variant.hgvs_p, Variant.gene_name, Variant.annotation, func.count().label('total'))\
            .join(Variant)
        if gene_name is not None:
            query = query.filter(and_(Variant.gene_name == gene_name, Variant.annotation != SYNONYMOUS_VARIANT))
        else:
            query = query.filter(Variant.annotation != SYNONYMOUS_VARIANT)
        query = query.group_by(
            VariantObservation.position, VariantObservation.reference,
            VariantObservation.alternate, Variant.hgvs_p, Variant.gene_name, Variant.annotation)
        query = query.order_by(desc('total')).limit(top)
        top_occurring_variants = pd.read_sql(query.statement, self.session.bind)

        if top_occurring_variants is not None and top_occurring_variants.shape[0] > 0:
            variant_counts_by_month = []
            # NOTE: one query per variant for the counts per month, will this be efficient?
            for _, variant in top_occurring_variants.iterrows():
                variant_counts_by_month.append(self.get_variant_counts_by_month(
                    variant.position, variant.reference, variant.alternate))
            top_occurring_variants_by_month = pd.concat(variant_counts_by_month)

            # join both tables with total counts and counts per month
            top_occurring_variants = top_occurring_variants.set_index(["position", "reference", "alternate"])\
                .join(
                top_occurring_variants_by_month.set_index(["position", "reference", "alternate"]))\
                .reset_index()

            # format the month column appropriately
            top_occurring_variants.month = top_occurring_variants.month.transform(
                lambda d: "{}-{:02d}".format(d.year, int(d.month)))

            # formats the DNA mutation
            top_occurring_variants['dna_mutation'] = top_occurring_variants.apply(
                lambda row: "{}:{}>{}".format(row["position"], row["reference"], row["alternate"]), axis=1)
            del top_occurring_variants["position"]
            del top_occurring_variants["reference"]
            del top_occurring_variants["alternate"]

            # replace the total count by the frequency
            count_samples = self.count_ena_samples_loaded()
            top_occurring_variants['frequency'] = top_occurring_variants.total.transform(
                lambda t: round(float(t) / count_samples, 3))
            del top_occurring_variants['total']

            # pivots the table over months
            top_occurring_variants = pd.pivot_table(
                top_occurring_variants, index=['gene_name', 'dna_mutation', 'hgvs_p', 'annotation', "frequency"],
                columns=["month"], values=["count"], fill_value=0).droplevel(0, axis=1).reset_index()

        return top_occurring_variants

    def get_variant_counts_by_month(self, position, reference, alternate) -> pd.DataFrame:
        query = self.session.query(
            VariantObservation.position, VariantObservation.reference,
            VariantObservation.alternate, func.date_trunc('month', SampleEna.first_created).label("month"),
            func.count().label("count"))\
            .join(Variant)\
            .join(SampleEna, VariantObservation.sample == SampleEna.run_accession)\
            .filter(and_(VariantObservation.position == position,
                         VariantObservation.reference == reference, VariantObservation.alternate == alternate,
                         Variant.annotation != SYNONYMOUS_VARIANT))\
            .group_by(VariantObservation.position, VariantObservation.reference,
                      VariantObservation.alternate, func.date_trunc('month', SampleEna.first_created))
        return pd.read_sql(query.statement, self.session.bind)

    def get_variants_cooccurrence_by_gene(self, gene_name, min_cooccurrence=5) -> pd.DataFrame:
        """
        Returns the full cooccurrence matrix of all non synonymous variants in a gene with at least
        min_occurrences occurrences.
        """
        variant_one = aliased(Variant)
        variant_two = aliased(Variant)
        query = self.session.query(VariantCooccurrence,
                                   variant_one.hgvs_p.label("hgvs_p_one"),
                                   variant_one.annotation.label("annotation_one"),
                                   variant_two.hgvs_p.label("hgvs_p_two"),
                                   variant_two.annotation.label("annotation_two")) \
            .filter(VariantCooccurrence.count >= min_cooccurrence)\
            .join(variant_one, and_(VariantCooccurrence.chromosome_one == variant_one.chromosome,
                                    VariantCooccurrence.position_one == variant_one.position,
                                    VariantCooccurrence.reference_one == variant_one.reference,
                                    VariantCooccurrence.alternate_one == variant_one.alternate)) \
            .join(variant_two, and_(VariantCooccurrence.chromosome_two == variant_two.chromosome,
                                    VariantCooccurrence.position_two == variant_two.position,
                                    VariantCooccurrence.reference_two == variant_two.reference,
                                    VariantCooccurrence.alternate_two == variant_two.alternate)) \
            .filter(and_(variant_one.gene_name == gene_name,
                         variant_one.annotation != SYNONYMOUS_VARIANT,
                         variant_two.gene_name == gene_name,
                         variant_two.annotation != SYNONYMOUS_VARIANT))
        data = pd.read_sql(query.statement, self.session.bind)

        def get_variant_id(x):
            return "{position}:{reference}>{alternate}".format(position=x[0], reference=x[1], alternate=x[2])

        full_matrix_with_annotations = None
        if data.shape[0] > 0:
            # build a single column identifying each variant
            data["variant_one"] = data[['position_one', 'reference_one', 'alternate_one']].apply(get_variant_id, axis=1)
            data["variant_two"] = data[['position_two', 'reference_two', 'alternate_two']].apply(get_variant_id, axis=1)

            # save annotations to a different dataframe
            annotations = data[["variant_one", "variant_two", "hgvs_p_one", "annotation_one", "hgvs_p_two",
                                "annotation_two"]]

            # delete other columns
            del data["chromosome_one"]
            del data["position_one"]
            del data["reference_one"]
            del data["alternate_one"]
            del data["chromosome_two"]
            del data["position_two"]
            del data["reference_two"]
            del data["alternate_two"]
            del data["hgvs_p_one"]
            del data["annotation_one"]
            del data["hgvs_p_two"]
            del data["annotation_two"]

            # from the sparse matrix builds in memory the complete matrix
            # TODO: would plotly improve its heatmap implementation so we don't need this memory misuse??
            all_variants = pd.concat([data.variant_one, data.variant_two], axis=0)
            all_variants = all_variants.sort_values().unique()
            empty_table = pd.DataFrame(
                index=pd.MultiIndex.from_product([all_variants, all_variants], names=["variant_one", "variant_two"]))
            empty_table["count"] = 0

            # adds values into empty table
            full_matrix = empty_table + data.set_index(["variant_one", "variant_two"])
            #full_matrix.fillna(0, inplace=True)
            full_matrix_with_annotations = full_matrix.join(annotations.set_index(["variant_one", "variant_two"]))
            full_matrix_with_annotations.reset_index(inplace=True)

            # removes the upper triangle
            full_matrix_with_annotations["count"] = full_matrix_with_annotations[[
                "variant_one", "variant_two", "count"]].apply(
                lambda x: x["count"] if x.variant_one <= x.variant_two else None, axis=1)

        return full_matrix_with_annotations
