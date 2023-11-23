from logzero import logger
from pathlib import Path
import pandas as pd
import numpy as np
from covigator.precomputations import GENES_DICT, POLYPROTEIN_OFFSET
from sklearn.linear_model import LinearRegression
import joblib
from sqlalchemy.orm import Session
from covigator.database.queries import Queries
from covigator.database.model import SubclonalVariantObservation


class FeatureLoader:
    DATA = Path(__file__) / "model_data"

    @staticmethod
    def _load_sift_features(intrahost_observations: pd.DataFrame) -> pd.DataFrame:
        """
        Load SIFT pre-computed annotation scores and annotate intrahost_observations
        """
        sift = pd.read_csv(FeatureLoader.DATA / 'Sars_cov_2.all_variants_SIFTannotations.tsv',
                           sep='\t',
                           usecols=['GENE_NAME', 'AMINO_POS', 'REF_AMINO', 'ALT_AMINO', 'SIFT_SCORE'])

        columns = ['variant_id', 'gene_name', 'position_amino_acid', 'reference_amino_acid', 'alternate_amino_acid']
        assert all([x in columns for x in intrahost_observations.columns]), "Missing columns in intrahost table"
        #intrahost_observations = pd.read_parquet(
        #    os.path.join(get_covigator_data_folder(), "subclonal_variant_observation_samples.parquet"),
        #    columns=columns)
        intrahost_observations_sift = intrahost_observations.join(
            sift.set_index(['GENE_NAME', 'AMINO_POS', 'REF_AMINO', 'ALT_AMINO']), \
            on=['gene_name', 'position_amino_acid', 'reference_amino_acid', 'alternate_amino_acid'])
        intrahost_observations_sift.rename(columns={'SIFT_SCORE': 'sift_score'}, inplace=True)

        sift_feature = intrahost_observations_sift.groupby('variant_id').first()
        sift_feature = sift_feature[~sift_feature.position_amino_acid.isna()]

        logger.info("Number of variants {} of which {} have a missing SIFT score".format(
            sift_feature.shape[0], sift_feature[sift_feature["sift_score"].isna()].shape[0]))

        logger.info(sift_feature.sift_score.describe())

        return sift_feature

    @staticmethod
    def _load_mutfunc_features(intrahost_observations: pd.DataFrame) -> pd.DataFrame:
        """
        Load mutfunc resoource and annotate intrahost_observations
        """
        mutfunc = pd.read_csv(FeatureLoader.DATA / 'mutfunc_sars/summary.tsv', sep='\t')

        # translate gene names
        mutfunc['gene_name'] = mutfunc['name'].transform(lambda x: GENES_DICT.get(x))

        mutfunc['gene_position'] = mutfunc[['name', 'position']].apply(
            lambda x: POLYPROTEIN_OFFSET.get(x[0], 0) + x[1], axis=1)

        # maps the functional annotations in protein coordinates to our variant ids
        # WARNING: we lose the annotations from all ORF1ab as these are annotated with protein domain coordinates...
        # who would have thought?
        columns = ['variant_id', 'gene_name', 'position_amino_acid', 'reference_amino_acid', 'alternate_amino_acid']
        assert all([x in columns for x in intrahost_observations.columns]), "Missing columns in intrahost table"
        #intrahost_observations = pd.read_parquet(
        #    os.path.join(get_covigator_data_folder(), "subclonal_variant_observation_samples.parquet"),
        #    columns=columns)
        intrahost_observations_mutfunc = intrahost_observations.join(
            mutfunc.set_index(['gene_name', 'gene_position', 'wt', 'mut']), \
            on=['gene_name', 'position_amino_acid', 'reference_amino_acid', 'alternate_amino_acid'])

        mutfunc_features = intrahost_observations_mutfunc.groupby('variant_id').first()
        mutfunc_features = mutfunc_features[~mutfunc_features.position_amino_acid.isna()]

        logger.info("Number of variants {}".format(mutfunc_features.shape[0]))
        logger.info(mutfunc_features.describe())

        return mutfunc_features

    @staticmethod
    def _build_features_table(
            variants: pd.DataFrame,
            features_vaf_regression: pd.DataFrame,
            features_conservation: pd.DataFrame,
            features_mutfunc: pd.DataFrame,
            feature_sift: pd.DataFrame,
            features_dnds: pd.DataFrame,
            feature_iedb: pd.DataFrame = None,
            fill_na=None):
        """
        variants is a DataFrame containing a column "variant_id"
        """

        features_table = variants[['variant_id']].set_index('variant_id') \
            .join(features_vaf_regression, how='left') \
            .join(features_conservation, how='left') \
            .join(features_mutfunc[['mut_escape_mean', 'mut_escape_max']], how='left') \
            .join(feature_sift[['sift_score']], how='left')

        if feature_iedb is not None:
            # NOTE: we decided to exclude IEDB features in most models
            features_table = features_table \
                .join(feature_iedb, how='left')

        features_table = features_table \
            .join(features_dnds[['clonal_dnds', 'subclonal_dnds']], how='left')

        # sorts table by variant_id
        features_table.sort_index(inplace=True)

        if fill_na is not None:
            # fill NA values with whatever value was provided, this is necessary at least for RF models
            features_table.fillna(fill_na, inplace=True)

        return features_table

    @staticmethod
    def _linear_regression(x, y, model):
        # TODO: do we want to add the intersection as an additional feature?
        #  that would represent how high is the starting point in terms of VAFs
        try:
            fit = model.fit(np.array(x).reshape((-1, 1)), np.array(y))
            return fit.coef_[0]
        except ValueError:
            return None
        except TypeError:
            return None

    @staticmethod
    def _calculate_vaf_regression(intrahost_observations_samples: pd.DataFrame) -> pd.DataFrame:
        """
        Input DF requires the columns "variant_id", "vaf" and "collection_date"
        Output DF is grouped by "variant_id" and has an additional column "vaf_regression_coefficient"
        """
        model = LinearRegression()
        intrahost_observations_samples['collection_date_timestamp'] = \
            intrahost_observations_samples[
                    ~intrahost_observations_samples.collection_date.isna()]['collection_date'].apply(
                        lambda x: x.timestamp())
        intrahost_regression_vafs = \
            intrahost_observations_samples[['variant_id', 'vaf', 'collection_date_timestamp']] \
                .groupby(['variant_id']) \
                .agg({'vaf': list, 'collection_date_timestamp': list})

        intrahost_regression_vafs['vaf_regression_coefficient'] = intrahost_regression_vafs.apply(
            lambda x: FeatureLoader._linear_regression(x[1], x[0], model), axis=1)

        return intrahost_regression_vafs[["vaf_regression_coefficient"]]

    @staticmethod
    def _extract_conservation_features(intrahost_observations_samples):
        """
        Input DF requires the columns 'variant_id', 'cons_hmm_sars_cov_2', 'cons_hmm_sarbecovirus'
        and 'cons_hmm_vertebrate_cov'
        Output DF is grouped by "variant_id" and contains the above input columns
        """
        conservation_features = intrahost_observations_samples[
            ['variant_id', 'cons_hmm_sars_cov_2', 'cons_hmm_sarbecovirus', 'cons_hmm_vertebrate_cov']] \
            .groupby('variant_id').first()
        return conservation_features


class IntrahostPredictionLoader:
    DATA = Path(__file__) / "model_data"

    def __init__(self, session: Session):
        self.model_file = IntrahostPredictionLoader.DATA / "05_model7_xgboost.joblib"
        self.model = joblib.load(self.model_file)
        self.session = session
        self.queries = Queries(session=self.session)

    def get_intrahost_mutations(self):
        """
        Implement query to database to obtain intrahost mutations
            - Merge with ENA sample metadata, required?
        Calculate dn/ds stuff
            - Requires porting Rangas R-code to python function
        """
        intrahost_mutations = self.session.query(SubclonalVariantObservation.variant_id,
                                                 SubclonalVariantObservation.gene_name,
                                                 SubclonalVariantObservation.variant_type,
                                                 SubclonalVariantObservation.vaf,
                                                 SubclonalVariantObservation.cons_hmm_sars_cov_2,
                                                 SubclonalVariantObservation.cons_hmm_sarbecovirus,
                                                 SubclonalVariantObservation.cons_hmm_vertebrate_cov)
        intrahost_mutations = pd.read_sql(intrahost_mutations.statement, self.session.bind)
        intrahost_mutations = intrahost_mutations.drop_duplicates()
        return intrahost_mutations

    def get_clonal_mutations(self):
        pass

    def calculate_dn_ds_across_all_time(self):
        """
        THis will replace the R code used to calculate dnds ratios
        """
        codons = pd.read_csv(self.DATA / "codons.csv")
        trinucleotide_ns_s = pd.read_csv(self.DATA / "trinucelotide_gene_NS_S.csv")
        pass
    def _build_feature_table(self):
        pass

