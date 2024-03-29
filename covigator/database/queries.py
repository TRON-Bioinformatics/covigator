from datetime import date, datetime
from typing import List, Union, Tuple

import pandas as pd
from logzero import logger
import sqlalchemy
from sqlalchemy import and_, desc, asc, func, String, DateTime, text
from sqlalchemy.engine.default import DefaultDialect
from sqlalchemy.orm import Session, aliased
from sqlalchemy.sql.sqltypes import NullType
from covigator import SYNONYMOUS_VARIANT
from covigator.database.model import DataSource, SampleEna, JobStatus, \
    VariantObservation, Gene, Variant, VariantCooccurrence, Conservation, \
    SubclonalVariantObservation, PrecomputedVariantsPerSample, PrecomputedSubstitutionsCounts, PrecomputedIndelLength, \
    VariantType, PrecomputedAnnotation, PrecomputedOccurrence, PrecomputedTableCounts, \
    PrecomputedVariantAbundanceHistogram, PrecomputedSynonymousNonSynonymousCounts, RegionType, Domain, \
    LastUpdate, SampleCovid19Portal, VariantObservationCovid19Portal, VariantCovid19Portal, Lineages, LineageVariant, \
    LineageDefiningVariants, NewsSection
from covigator.exceptions import CovigatorQueryException, CovigatorDashboardMissingPrecomputedData


class Queries:

    def __init__(self, session: Session):
        self.session = session

    @staticmethod
    def get_variant_observation_klass(source: str):
        if source == DataSource.ENA.name:
            klass = VariantObservation
        elif source == DataSource.COVID19_PORTAL.name:
            klass = VariantObservationCovid19Portal
        else:
            raise CovigatorQueryException("Bad data source: {}".format(source))
        return klass

    @staticmethod
    def get_variant_klass(source: str):
        if source == DataSource.ENA.name:
            klass = Variant
        elif source == DataSource.COVID19_PORTAL.name:
            klass = VariantCovid19Portal
        else:
            raise CovigatorQueryException("Bad data source: {}".format(source))
        return klass

    @staticmethod
    def get_sample_klass(source: str):
        if source == DataSource.ENA.name:
            klass = SampleEna
        elif source == DataSource.COVID19_PORTAL.name:
            klass = SampleCovid19Portal
        else:
            raise CovigatorQueryException("Bad data source: {}".format(source))
        return klass

    @staticmethod
    def get_variant_cooccurrence_klass(source: str):
        if source == DataSource.ENA.name:
            klass = VariantCooccurrence
        else:
            raise CovigatorQueryException("Bad data source: {}".format(source))
        return klass

    def find_job_by_accession_and_status(
            self, run_accession: str, status: JobStatus, data_source: DataSource) -> Union[SampleEna]:
        klass = self.get_sample_klass(source=data_source.name)
        return self.session.query(klass)\
            .filter(and_(klass.run_accession == run_accession, klass.status == status)).first()

    def find_job_by_accession(self, run_accession: str, data_source: DataSource) -> Union[SampleEna]:
        klass = self.get_sample_klass(source=data_source.name)
        return self.session.query(klass).filter(klass.run_accession == run_accession).first()

    def find_first_by_status(self, data_source: DataSource, status: Tuple, n=100) -> List[Union[SampleEna]]:
        klass = self.get_sample_klass(source=data_source.name)
        return self.session.query(klass) \
            .filter(klass.status.in_(status)) \
            .order_by(klass.created_at.desc()) \
            .limit(n) \
            .all()

    def find_first_pending_jobs(
            self, data_source: DataSource, n=100, status: Tuple = (JobStatus.DOWNLOADED, )) -> List[Union[SampleEna]]:
        return self.find_first_by_status(data_source=data_source, status=status, n=n)

    def find_first_jobs_to_download(self, data_source: DataSource, n=100) -> List[Union[SampleEna]]:
        return self.find_first_by_status(data_source=data_source, status=(JobStatus.PENDING, ), n=n)

    def count_jobs_in_queue(self, data_source):
        return self.count_jobs_by_status(data_source=data_source, status=JobStatus.QUEUED)

    def count_jobs_by_status(self, data_source: DataSource, status: JobStatus):
        klass = self.get_sample_klass(source=data_source.name)
        return self.session.query(klass).filter(klass.status == status).count()

    def find_sample_by_accession(self, run_accession: str, source: DataSource) -> Union[SampleEna]:
        klass = self.get_sample_klass(source=source.name)
        return self.session.query(klass).filter(klass.run_accession == run_accession).first()

    def get_countries(self, source: str) -> List[str]:
        klass = self.get_sample_klass(source=source)
        return [c for c, in self.session.query(klass.country).filter(
                klass.status == JobStatus.FINISHED.name).distinct().order_by(klass.country.asc()).all()]

    def get_lineages(self, source: str) -> List[str]:
        klass = self.get_sample_klass(source=source)
        return [c for c, in self.session.query(klass.pangolin_lineage).filter(
            and_(klass.status == JobStatus.FINISHED.name, klass.pangolin_lineage != None)).distinct().order_by(
                klass.pangolin_lineage.asc()).all()]

    def find_parent_who_label(self, lineage, parent_mapping) -> str:
        """
        In some instances a constellation does not include a WHO label. Nevertheless, these sublineages also belong to
        this VOC and should also be grouped with this label. Returns a string containing the WHO label of the topmost
        parent in the tree if present, otherwise just the label.

        Examples are local lineage variations such as AY.4.2 and AY.4 --> Delta
        """
        lineage = parent_mapping.query("pangolin_lineage == @lineage").get(["pangolin_lineage", "who_label",
                                                                            "parent_lineage_id"])
        # No or last parent in tree --> return WHO label
        parent = lineage.get("parent_lineage_id").item()
        if pd.isnull(parent):
            return lineage.get("who_label").item()
        else:
            return self.find_parent_who_label(parent, parent_mapping)

    def get_lineages_who_label(self) -> pd.DataFrame:
        """
        Query database for lineage WHO label annotation. Returns a DataFrame with columns: pangolin_lineage, who_label
        """
        query = self.session.query(Lineages.pango_lineage_id.label("pangolin_lineage"), Lineages.who_label)
        lineages = pd.read_sql(query.statement, self.session.bind)
        lineages = lineages[["pangolin_lineage", "who_label"]]
        return lineages

    def get_combined_labels(self, source: str) -> pd.DataFrame:
        """
        Create a mapping from pangolin IDs in the data source to WHO identifiers and create a combined
        label to be used in the dashboard
        """
        who_labels = self.get_lineages_who_label()
        lineages = self.get_lineages(source)
        if not lineages:
            return None
        lineages = pd.DataFrame(lineages, columns=["pangolin_lineage"])
        lineages = lineages.merge(who_labels, how="left")
        lineages["combined_label"] = lineages.apply(lambda x: f"{x.pangolin_lineage} - {x.who_label}"
            if not pd.isnull(x.who_label) else f"{x.pangolin_lineage}", axis=1)
        return lineages

    def get_lineage_defining_variants(self) -> pd.DataFrame:
        """
        Query database for lineage defining variants. Returns a dataframe with columns: pangolin_lineage, variant
        """
        query = self.session.query(LineageVariant.pango_lineage_id.label("pangolin_lineage"), LineageVariant.variant_id,
                                   LineageDefiningVariants.hgvs).join(LineageDefiningVariants)
        lineage_variants = pd.read_sql(query.statement, self.session.bind)
        # Rename columns to match columns used in recurrent/intrahost mutations tab
        lineage_variants["dna_mutation"] = lineage_variants.apply(lambda x: x.variant_id if pd.isnull(x.hgvs)
                                                                  else None, axis=1)
        lineage_variants = lineage_variants.rename(columns={"hgvs": "hgvs_p"})
        aa_level_mutations = lineage_variants[~pd.isnull(lineage_variants.hgvs_p)].loc[:, ['pangolin_lineage', 'hgvs_p']]
        aa_level_mutations = aa_level_mutations.sort_values(['hgvs_p', 'pangolin_lineage']).groupby('hgvs_p')['pangolin_lineage'].agg(','.join)
        aa_level_mutations = aa_level_mutations.reset_index()

        nucleotide_level_mutations = lineage_variants[~pd.isnull(lineage_variants.dna_mutation)].loc[:, ['pangolin_lineage','dna_mutation']]
        nucleotide_level_mutations = nucleotide_level_mutations.sort_values(['dna_mutation', 'pangolin_lineage']).groupby('dna_mutation')['pangolin_lineage'].agg(','.join)
        nucleotide_level_mutations = nucleotide_level_mutations.reset_index()

        return aa_level_mutations, nucleotide_level_mutations

    def _merge_with_lineage_defining_variants(self, data: pd.DataFrame):
        """
        Merge tables from recurrent and intrahost mutations tab with lineage defining mutations
        """

        assert "variant_id" in data.columns, "Column variant_id is missing..."
        assert "hgvs_p" in data.columns, "Column hgvs_p is missing..."
        lineage_mutation_aa, lineage_mutation_nuc = self.get_lineage_defining_variants()
        lineage_mutation_nuc.rename(columns={'dna_mutation': 'variant_id'}, inplace=True)

        # Merge data with lineage defining variants on different levels
        data = data.merge(lineage_mutation_aa, how="left", left_on="hgvs_p", right_on="hgvs_p")
        data = data.merge(lineage_mutation_nuc, how="left", left_on="variant_id", right_on="variant_id")

        data["pangolin_lineage_x"].fillna(data["pangolin_lineage_y"], inplace=True)
        data.drop(columns=["pangolin_lineage_y"], inplace=True)
        # formats the lineage column
        data.rename(columns={"pangolin_lineage_x": "pangolin_lineage"}, inplace=True)
        data["no_of_lineages"] = data.fillna({"pangolin_lineage": ""}).apply(
                lambda x: len(x.pangolin_lineage.split(",")), axis=1)
        data['pangolin_hover'] = data.apply(
                lambda x: "{} lineages".format(x.no_of_lineages) if x.no_of_lineages > 3 else x.pangolin_lineage,
                    axis=1)
        return data

    def get_variants_per_sample(self, data_source: str, genes: List[str], variant_types: List[str]):
        """
        Returns a DataFrame with columns: number_mutations, count, type
        where type: SNV, insertion or deletion
        """
        self._assert_data_source(data_source)
        query = self.session.query(PrecomputedVariantsPerSample)\
            .filter(PrecomputedVariantsPerSample.source == data_source)
        if genes is not None and genes:
            query = query.filter(PrecomputedVariantsPerSample.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedVariantsPerSample.gene_name == None)
        if variant_types is not None and variant_types:
            query = query.filter(PrecomputedVariantsPerSample.variant_type.in_(variant_types))

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            data.variant_type = data.variant_type.transform(lambda x: x.name if x else x)
        return data[["number_mutations", "variant_type", "count"]] \
            .groupby(["number_mutations", "variant_type"]).sum().reset_index()

    def _assert_data_source(self, data_source):
        if data_source != DataSource.ENA.name and data_source != DataSource.COVID19_PORTAL.name:
            raise CovigatorQueryException("Bad data source:  {}".format(data_source))

    def get_indel_lengths(self, data_source, genes):
        self._assert_data_source(data_source)
        query = self.session.query(PrecomputedIndelLength).filter(PrecomputedIndelLength.source == data_source)
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
        self._assert_data_source(data_source)
        query = self.session.query(PrecomputedAnnotation).filter(PrecomputedAnnotation.source == data_source)
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
        self._assert_data_source(data_source)
        query = self.session.query(PrecomputedSubstitutionsCounts)\
            .filter(PrecomputedSubstitutionsCounts.source == data_source)
        if variant_types is not None and variant_types:
            query = query.filter(PrecomputedSubstitutionsCounts.variant_type.in_(variant_types))
        if genes is not None and genes:
            query = query.filter(PrecomputedSubstitutionsCounts.gene_name.in_(genes))
        else:
            query = query.filter(PrecomputedSubstitutionsCounts.gene_name == None)

        data = pd.read_sql(query.statement, self.session.bind)
        if data.shape[0] > 0:
            data.variant_type = data.variant_type.transform(lambda x: x.name if x else x)
            data["substitution"] = data[["reference", "alternate"]].apply(lambda x: "{}>{}".format(x[0], x[1]), axis=1)
            total_number_variants = data["count"].sum()
            data = data[["substitution", "variant_type", "count"]] \
                .groupby(["substitution", "variant_type"]).sum().reset_index()\
                .sort_values("count", ascending=False)
            data["rate"] = (data["count"] / total_number_variants).transform(lambda x: "{} %".format(round(x * 100, 1)))
            data = data.head(12)
        return data

    def get_accumulated_samples_by_country(
            self, data_source: str, countries: List[str], lineages: List[str], min_samples=100) -> pd.DataFrame:
        """
        Returns a DataFrame with columns: data, country, cumsum, count
        """
        klass = self.get_sample_klass(source=data_source)
        query = self.session.query(
            func.count().label("count"), klass.collection_date.label("date"), klass.country) \
            .filter(klass.status == JobStatus.FINISHED.name) \
            .group_by(klass.collection_date, klass.country)
        if countries:
            query = query.filter(klass.country.in_(countries))
        if lineages:
            query = query.filter(klass.pangolin_lineage.in_(lineages))
        samples = pd.read_sql(query.statement, self.session.bind).astype({'date': 'datetime64', 'count': 'float64'})

        filled_table = None
        if samples is not None and samples.shape[0] > 0:
            # accumulates count ordered by date
            samples['cumsum'] = samples.groupby(['country'])['count'].cumsum()

            # creates empty table with all pairwise combinations of date and country
            dates = samples.date[~samples.date.isna()].sort_values().unique()
            countries = samples[(~samples.country.isna()) & (samples["cumsum"] >= min_samples)] \
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

    def get_accumulated_lineages_by_country(
            self, data_source: str, countries: List[str], lineages: List[str], min_samples=100) -> pd.DataFrame:
        """
        Returns a DataFrame with columns: data, country, cumsum, count
        """
        klass = self.get_sample_klass(source=data_source)
        query = self.session.query(
            func.count().label("count"), klass.collection_date.label("date"), klass.pangolin_lineage.label("lineage")) \
            .filter(and_(klass.status == JobStatus.FINISHED.name, klass.pangolin_lineage != None)) \
            .group_by(klass.collection_date, klass.pangolin_lineage)
        if countries:
            query = query.filter(klass.country.in_(countries))
        if lineages:
            query = query.filter(klass.pangolin_lineage.in_(lineages))
        samples = pd.read_sql(query.statement, self.session.bind).astype({'date': 'datetime64', 'count': 'float64'})

        filled_table = None
        if samples is not None and samples.shape[0] > 0:
            # accumulates count ordered by date
            samples['cumsum'] = samples.groupby(['lineage'])['count'].cumsum()

            # creates empty table with all pairwise combinations of date and country
            dates = samples.date[~samples.date.isna()].sort_values().unique()

            lineages = samples[(~samples.lineage.isna()) & (samples["cumsum"] >= min_samples)] \
                .sort_values("cumsum", ascending=False).lineage.unique()

            empty_table = pd.DataFrame(
                index=pd.MultiIndex.from_product([dates, lineages], names=["date", "lineage"]))
            empty_table["count"] = 0

            # adds values into empty table
            filled_table = pd.merge(
                left=empty_table.reset_index(),
                right=samples.loc[:, ["date", "lineage", "count"]],
                on=["date", "lineage"],
                how='left',
                suffixes=("_x", "_y")).fillna(0)

            filled_table["count"] = filled_table.count_x + filled_table.count_y
            filled_table['cumsum'] = filled_table.groupby(['lineage'])['count'].cumsum()

            # add total samples per day
            counts_per_date = filled_table[["date", "count"]].groupby("date").sum().reset_index()
            counts_per_date.rename(columns={'count': 'total_per_date'}, inplace=True)
            filled_table = pd.merge(
                left=filled_table,
                right=counts_per_date,
                on="date")
            filled_table['ratio_per_date'] = filled_table[["count", "total_per_date"]].apply(
                lambda x: float(x[0]) / x[1], axis=1)

        return filled_table

    def get_sample_months(self, pattern, data_source: str) -> List[datetime]:
        klass = self.get_sample_klass(source=data_source)
        dates = [d.strftime(pattern) for d, in self.session.query(klass.collection_date).filter(
            and_(klass.status == JobStatus.FINISHED.name, klass.collection_date.isnot(None))).distinct().all()]
        return sorted(list(set(dates)))

    def get_gene(self, gene_name: str) -> Gene:
        return self.session.query(Gene).filter(Gene.name == gene_name).first()

    
    def get_genes(self) -> List[Gene]:
        return self.session.query(Gene).order_by(Gene.start).all()

    
    def get_genes_df(self) -> pd.DataFrame:
        return pd.read_sql(self.session.query(Gene).order_by(Gene.start).statement, self.session.bind)

    
    def get_domains_df(self) -> pd.DataFrame:
        return pd.read_sql(self.session.query(Domain).order_by(Domain.start).statement, self.session.bind)

    
    def get_domain(self, domain_name: str) -> Domain:
        return self.session.query(Domain).filter(Domain.name == domain_name).first()

    
    def get_domains(self) -> List[Domain]:
        return self.session.query(Domain).order_by(Domain.gene_name, Domain.start).all()

    
    def get_domains_by_gene(self, gene_name: str) -> List[Domain]:
        return self.session.query(Domain).filter(Domain.gene_name == gene_name).order_by(Domain.start).all()

    
    def get_non_synonymous_variants_by_region(self, start, end, source) -> pd.DataFrame:

        klass = self.get_variant_observation_klass(source)
        query = self.session.query(
            klass.position, klass.annotation_highest_impact, klass.hgvs_p, func.count().label("count_occurrences"))\
            .filter(and_(klass.annotation_highest_impact != SYNONYMOUS_VARIANT,
                         klass.position >= start, klass.position <= end))

        subquery = query.group_by(klass.position, klass.annotation_highest_impact, klass.hgvs_p).subquery()
        return pd.read_sql(
            self.session.query(subquery).filter(subquery.c.count_occurrences > 1).statement, self.session.bind)
    
    def count_samples(self, source: str, cache=True) -> int:
        self._assert_data_source(source)
        if cache:
            query = self.session.query(PrecomputedTableCounts.count)
            if source == DataSource.ENA.name:
                query = query.filter(and_(
                    PrecomputedTableCounts.table == SampleEna.__name__,
                    PrecomputedTableCounts.factor == PrecomputedTableCounts.FACTOR_SOURCE
                ))
            elif source == DataSource.COVID19_PORTAL.name:
                query = query.filter(and_(
                    PrecomputedTableCounts.table == SampleCovid19Portal.__name__,
                    PrecomputedTableCounts.factor == PrecomputedTableCounts.FACTOR_SOURCE
                ))
            else:
                raise ValueError("Unknown data source")
            result = query.first()
            if result is None:
                raise CovigatorDashboardMissingPrecomputedData
            count = result.count
        else:
            klass = self.get_sample_klass(source=source)
            count = self.session.query(klass).filter(klass.status == JobStatus.FINISHED.name).count()
        return count

    def count_countries(self, source: str = None, cache=True):
        if cache:
            query = self.session.query(PrecomputedTableCounts.count).filter(and_(
                    PrecomputedTableCounts.table == PrecomputedTableCounts.VIRTUAL_TABLE_COUNTRY,
                    PrecomputedTableCounts.factor == PrecomputedTableCounts.FACTOR_SOURCE,
                    PrecomputedTableCounts.value == source
                ))
            result = query.first()
            if result is None:
                raise CovigatorDashboardMissingPrecomputedData
            count = result.count
        else:
            count = len(self.get_countries(source=source))
        return count

    def count_variants(self, source: str, cache=True):
        klass = self.get_variant_klass(source=source)
        if cache:
            result = self.session.query(PrecomputedTableCounts.count) \
                .filter(and_(
                PrecomputedTableCounts.table == klass.__name__, PrecomputedTableCounts.value == source)).first()
            if result is None:
                raise CovigatorDashboardMissingPrecomputedData
            count = result.count
        else:
            count = self.session.query(klass).count()
        return count

    def count_variant_observations(self, source: str = None, cache=True):
        klass = self.get_variant_observation_klass(source)
        if cache:
            result = self.session.query(PrecomputedTableCounts.count).filter(
                PrecomputedTableCounts.table == klass.__name__).first()
            if result is None:
                raise CovigatorDashboardMissingPrecomputedData
            count = result.count
        else:
            count = self.session.query(klass).count()
        return count

    def count_subclonal_variant_observations(self, cache=True):
        if cache:
            query = self.session.query(PrecomputedTableCounts.count) \
                .filter(PrecomputedTableCounts.table == SubclonalVariantObservation.__name__)
            result = query.first()
            count = result.count
        else:
            count = self.session.query(SubclonalVariantObservation).count()
        return count

    def count_unique_subclonal_variant(self, cache=True):
        if cache:
            query = self.session.query(PrecomputedTableCounts.count) \
                .filter(PrecomputedTableCounts.table == SubclonalVariantObservation.__name__ + "_unique")
            result = query.first()
            count = result.count
        else:
            count = self.session.query(SubclonalVariantObservation)\
                .distinct(SubclonalVariantObservation.variant_id).count()
        return count

    def count_unique_only_subclonal_variant(self, cache=True):
        if cache:
            query = self.session.query(PrecomputedTableCounts.count) \
                .filter(PrecomputedTableCounts.table == SubclonalVariantObservation.__name__ + "_unique_only_subclonal")
            result = query.first()
            count = result.count
        else:
            sql_query = """
                    SELECT COUNT(*) as count from (
                        select distinct variant_id 
                        from {subclonal_variants_table}
                        where variant_id not in (
                            select distinct variant_id from {variant_observations_table}
                        )
                        and vaf >= 0.03
                    ) as variants""".format(
                subclonal_variants_table=SubclonalVariantObservation.__tablename__,
                variant_observations_table=VariantObservation.__tablename__)
            count = int(pd.read_sql_query(sql_query, self.session.bind)["count"][0])
        return count

    def get_date_of_first_sample(self, source: DataSource = DataSource.ENA) -> date:
        """
        Returns the date of the earliest ENA sample loaded in the database
        """
        klass = self.get_sample_klass(source=source.name)
        result = self.session.query(klass.collection_date).filter(
            and_(klass.status == JobStatus.FINISHED.name, klass.collection_date.isnot(None))).order_by(asc(klass.collection_date)).first()
        return result[0] if result is not None else result

    def get_date_of_most_recent_sample(self, source: DataSource = DataSource.ENA) -> date:
        """
        Returns the date of the latest ENA sample loaded in the database
        """
        klass = self.get_sample_klass(source=source.name)
        result = self.session.query(klass.collection_date).filter(
            and_(klass.status == JobStatus.FINISHED.name, klass.collection_date.isnot(None))).order_by(desc(klass.collection_date)).first()
        return result[0] if result is not None else result

    def get_last_update(self, data_source: DataSource) -> date:
        result = self.session.query(LastUpdate.update_time).filter(LastUpdate.source == data_source).order_by(
            desc(LastUpdate.update_time)).first()
        most_recent_update = result[0] if result is not None else result
        return most_recent_update

    def get_variant_counts_by_month(self, variant_id, source: str) -> pd.DataFrame:

        klass = self.get_variant_observation_klass(source=source)
        sample_klass = self.get_sample_klass(source=source)
        sql_query_ds_ena = text("""
        select count(*) as count, variant_id, date_trunc('month', date::timestamp) as month 
            from {variant_observation_table} 
            where variant_id=:variant_id
            and sample in (select run_accession from {sample_table} where status='FINISHED')
            group by variant_id, date_trunc('month', date::timestamp);
            """.format(
            variant_observation_table=klass.__tablename__,
            sample_table=sample_klass.__tablename__,
        ))
        data = pd.read_sql_query(sql_query_ds_ena, self.session.bind, params={"variant_id": variant_id})
        data['month'] = pd.to_datetime(data['month'], utc=True)
        return data[~data.month.isna()]

    def get_sample_counts_by_month(self, source: str) -> pd.DataFrame:
        klass = self.get_sample_klass(source=source)
        # NOTE: this query was originally implemented with SQLAlchemy syntax, but the func.date_trunc function
        # provides different results. Do not change back!
        query = """
        select date_trunc('month', collection_date::timestamp) as month,
            count(*) as sample_count
            from {table}
            where status='FINISHED'
            group by date_trunc('month', collection_date::timestamp);
            """.format(
            table=klass.__tablename__
        )
        counts = pd.read_sql(text(query), self.session.bind)
        counts['month'] = pd.to_datetime(counts['month'], utc=True)
        return counts

    def get_sparse_cooccurrence_matrix(self, gene_name, domain, source: str, min_cooccurrence=5) -> pd.DataFrame:
        """
        Returns the sparse cooccurrence matrix of all non synonymous variants in a gene or domain with at least
        min_occurrences occurrences.
        """
        # query for cooccurrence matrix
        variant_cooccurrence_klass = self.get_variant_cooccurrence_klass(source=source)
        variant_one = aliased(Variant)
        variant_two = aliased(Variant)
        query = self.session.query(variant_cooccurrence_klass,
                                   variant_one.position,
                                   variant_one.reference,
                                   variant_one.alternate,
                                   variant_one.hgvs_p,
                                   variant_one.hgvs_p.label("hgvs_p_one"),
                                   variant_two.hgvs_p.label("hgvs_p_two"),
                                   (variant_one.hgvs_p + " - " + variant_two.hgvs_p).label("hgvs_tooltip"))

        query = query.filter(variant_cooccurrence_klass.count >= min_cooccurrence) \
            .join(variant_one, and_(variant_cooccurrence_klass.variant_id_one == variant_one.variant_id)) \
            .join(variant_two, and_(variant_cooccurrence_klass.variant_id_two == variant_two.variant_id))

        if domain is not None:
            query = query.filter(and_(
                variant_one.pfam_name == domain,
                variant_two.pfam_name == domain))
        elif gene_name is not None:
            query = query.filter(and_(
                variant_one.gene_name == gene_name,
                variant_two.gene_name == gene_name))

        self._print_query(query=query)
        data = pd.read_sql(query.statement, self.session.bind)

        return data

    def get_top_occurring_variants_precomputed(
            self, top=10, gene_name=None, domain=None, metric="count", source=None) -> pd.DataFrame:
        """
        Returns the top occurring variants + the segregated counts of occurrences per month
        with columns: chromosome, position, reference, alternate, total, month, count
        """
        query = self.session.query(PrecomputedOccurrence).filter(PrecomputedOccurrence.source == source)
        # if domain is provided it supersedes the gene filter
        if domain is not None:
            query = query.filter(PrecomputedOccurrence.domain == domain)
        elif gene_name is not None:
            query = query.filter(PrecomputedOccurrence.gene_name == gene_name)
        if metric == "count":
            query = query.order_by(PrecomputedOccurrence.count.desc())
        elif metric == "frequency_by_month":
            query = query.order_by(PrecomputedOccurrence.frequency.desc())
        else:
            raise CovigatorQueryException("Not supported metric for top occurring variants")
        
        top_occurring_variants = pd.read_sql(query.statement, self.session.bind)
        # Merge with lineage defining variants
        top_occurring_variants = self._merge_with_lineage_defining_variants(top_occurring_variants)
        # formats the DNA mutation
        top_occurring_variants.rename(columns={"variant_id": "dna_mutation"}, inplace=True)

        # pivots the table over months
        top_occurring_variants = pd.pivot_table(
            top_occurring_variants, index=['gene_name', 'dna_mutation', 'hgvs_p', 'annotation', "frequency", "total", "pangolin_lineage", "pangolin_hover"],
            columns=["month"], values=[metric], fill_value=0).droplevel(0, axis=1).reset_index()

        return top_occurring_variants.sort_values(by="frequency", ascending=False).head(top)

    def get_variant_abundance_histogram(self, source: str, bin_size=50, cache=True) -> pd.DataFrame:
        histogram = None
        if cache:
            query = self.session.query(PrecomputedVariantAbundanceHistogram)
            if source is not None:
                query = query.filter(and_(PrecomputedVariantAbundanceHistogram.bin_size == bin_size,
                                          PrecomputedVariantAbundanceHistogram.source == source))
            else:
                query = query.filter(and_(PrecomputedVariantAbundanceHistogram.bin_size == bin_size,
                                          PrecomputedVariantAbundanceHistogram.source == None))
            histogram = pd.read_sql(query.statement, self.session.bind)
            if histogram.shape[0] == 0:
                raise CovigatorDashboardMissingPrecomputedData
            histogram = histogram[["position_bin", "count_unique_variants", "count_variant_observations"]]
        else:
            klass = self.get_variant_klass(source=source)
            klass_observation = self.get_variant_observation_klass(source=source)

            # queries for the maximum position
            maximum_position = self.session.query(func.max(klass.position)).first()[0]
            if maximum_position is not None:
                # builds all possible bins
                all_bins = pd.DataFrame(
                    data=[i*bin_size for i in range(int(maximum_position/bin_size) + 1)],
                    columns=["position_bin"])

                # counts variants over those bins
                sql_query = text("""
                        SELECT cast("position"/:bin_size as int)*:bin_size AS position_bin,
                               COUNT(*) as count_unique_variants
                        FROM {table_name}
                        GROUP BY position_bin
                        ORDER BY position_bin;
                        """.format(table_name=klass.__tablename__))
                binned_counts_variants = pd.read_sql_query(sql_query, self.session.bind, params={'bin_size': bin_size})

                # counts variant observations over those bins
                sql_query = text("""
                        SELECT cast("position"/:bin_size as int)*:bin_size AS position_bin,
                               COUNT(*) as count_variant_observations
                        FROM {table_name}
                        GROUP BY position_bin
                        ORDER BY position_bin;
                        """.format(table_name=klass_observation.__tablename__))
                binned_counts_variant_observations = pd.read_sql_query(sql_query, self.session.bind,
                    params={'bin_size': bin_size})

                histogram = pd.merge(
                    left=pd.merge(
                        left=all_bins,
                        right=binned_counts_variants,
                        on="position_bin", how="left").reset_index(),
                    right=binned_counts_variant_observations,
                    on="position_bin", how="left").reset_index()
                histogram.fillna(0, inplace=True)
                histogram = histogram[["position_bin", "count_unique_variants", "count_variant_observations"]]

        return histogram

    def get_conservation_table(self, bin_size=50, start=None, end=None) -> pd.DataFrame:
        # counts variants over those bins
        sql_query = text("""
                SELECT cast("start"/:bin_size as int)*:bin_size AS position_bin,
                       AVG("conservation") as conservation,
                       AVG("conservation_sarbecovirus") as conservation_sarbecovirus,
                       AVG("conservation_vertebrates") as conservation_vertebrates
                FROM {table_name}
                {where}
                GROUP BY position_bin
                ORDER BY position_bin;
                """.format(table_name= Conservation.__tablename__,
                           where="WHERE start >= :start and start <= :end"
                           if start is not None and end is not None else ""))
        return pd.read_sql_query(sql_query, self.session.bind,
            params={"bin_size": bin_size, "start": start, "end": end})

    def get_dnds_table(self, source: str, countries=None, genes=None) -> pd.DataFrame:
        self._assert_data_source(data_source=source)
        # counts variants over those bins
        query_genes = self.session.query(PrecomputedSynonymousNonSynonymousCounts)\
            .filter(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.GENE)
        query_domains = self.session.query(PrecomputedSynonymousNonSynonymousCounts) \
            .filter(PrecomputedSynonymousNonSynonymousCounts.region_type == RegionType.DOMAIN)
        data_domains = None

        if source is not None:
            query_genes = query_genes.filter(PrecomputedSynonymousNonSynonymousCounts.source == source)
            query_domains = query_domains.filter(PrecomputedSynonymousNonSynonymousCounts.source == source)
        if countries is not None and len(countries) > 0:
            query_genes = query_genes.filter(PrecomputedSynonymousNonSynonymousCounts.country.in_(countries))
            query_domains = query_domains.filter(PrecomputedSynonymousNonSynonymousCounts.country.in_(countries))
        if genes is not None and len(genes) > 0:
            query_genes = query_genes.filter(PrecomputedSynonymousNonSynonymousCounts.region_name.in_(genes))
            domains = []
            for g in genes:
                domains = domains + [d.name for d in self.get_domains_by_gene(gene_name=g)]
            query_domains = query_domains.filter(PrecomputedSynonymousNonSynonymousCounts.region_name.in_(domains))
            data_domains = pd.read_sql(query_domains.statement, self.session.bind)

        data = pd.read_sql(query_genes.statement, self.session.bind)
        if data_domains is not None and data_domains.shape[1] > 0:
            data = pd.concat([data, data_domains])

        return data

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

    def column_windows(self, session, column, windowsize):
        """Return a series of WHERE clauses against
        a given column that break it into windows.

        Result is an iterable of tuples, consisting of
        ((start, end), whereclause), where (start, end) are the ids.

        Requires a database that supports window functions,
        i.e. Postgresql, SQL Server, Oracle.

        Enhance this yourself !  Add a "where" argument
        so that windows of just a subset of rows can
        be computed.

        """
        def int_for_range(start_id, end_id):
            if end_id:
                return and_(
                    column >= start_id,
                    column < end_id
                )
            else:
                return column >= start_id

        q = session.query(
            column,
            func.row_number(). \
                over(order_by=column). \
                label('rownum')).from_self(column)
        if windowsize > 1:
            q = q.filter(sqlalchemy.text("rownum %% %d=1" % windowsize))

        intervals = [id for id, in q]

        while intervals:
            start = intervals.pop(0)
            if intervals:
                end = intervals[0]
            else:
                end = None
            yield int_for_range(start, end)

    def windowed_query(self, query, column, windowsize):
        """"
        Break a Query into windows on a given column.

        This magic comes from here: https://github.com/sqlalchemy/sqlalchemy/wiki/RangeQuery-and-WindowedRangeQuery
        """
        for whereclause in self.column_windows(query.session, column, windowsize):
            for row in query.filter(whereclause).order_by(column):
                yield row

    def get_top_news(self, n=3):
        """
        Get news from database, sort by publishing date and return n news when specified.
        """
        query = self.session.query(NewsSection).order_by(NewsSection.published_date.desc())
        news = pd.read_sql(query.statement, self.session.bind)
        if n is not None:
            news = news.head(n)
        return news
