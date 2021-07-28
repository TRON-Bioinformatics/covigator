import functools

import pandas as pd
import numpy as np
import dash_core_components as dcc
import dash_table
from sqlalchemy.orm import Session
from covigator.dashboard.figures.figures import Figures
from covigator.database.model import SubclonalVariantObservation, VariantObservation, SampleEna
from covigator.exceptions import CovigatorQueryException


class SubclonalVariantsQueries:

    def __init__(self, session: Session):
        self.session = session

    def get_top_occurring_subclonal_variants(self, top, gene_name, min_vaf, order_by):

        if order_by == "score":
            order_by_clause = "order by score asc"
        elif order_by == "count":
            order_by_clause = "order by count_observations desc"
        elif order_by == "conservation":
            order_by_clause = "order by conservation asc"
        elif order_by == "vaf":
            order_by_clause = "order by median_vaf desc"
        else:
            raise CovigatorQueryException("Not supported order by of subclonal variants")

        # counts variants over those bins
        sql_query = """
        select variant_id, date, hgvs_p, gene_name, annotation_highest_impact, 
            cons_hmm_sars_cov_2, cons_hmm_sarbecovirus, cons_hmm_vertebrate_cov, pfam_name, pfam_description, vaf
            from {subclonal_variants_table_name} 
            where 
                vaf >= {min_vaf} and 
                sample not in (select run_accession from {samples_table_name} where library_strategy = 'RNA-Seq') and
                variant_id in (
                    select t.variant_id from (
                        select 
                            count(*) as count_observations, 
                            max(cons_hmm_sars_cov_2) as conservation, 
                            (ln(count(*)) * max(cons_hmm_sars_cov_2))::float as score, 
                            percentile_disc(0.5) within group (order by vaf) as median_vaf,
                            variant_id
                        from {subclonal_variants_table_name} 
                        where 
                            vaf >= {min_vaf} and 
                            gene_name is not null and 
                            {where_gene}
                            annotation_highest_impact != 'synonymous_variant' and 
                            sample not in (select run_accession from {samples_table_name} where library_strategy = 'RNA-Seq') 
                            and variant_id not in (
                                select distinct variant_id 
                                from {variants_table_name} 
                                where source='ENA' and gene_name is not null and annotation_highest_impact != 'synonymous_variant'
                            )
                        group by variant_id
                        {order_by}
                        limit {top}) as t);
            """.format(
            min_vaf=min_vaf,
            subclonal_variants_table_name=SubclonalVariantObservation.__tablename__,
            variants_table_name=VariantObservation.__tablename__,
            samples_table_name=SampleEna.__tablename__,
            top=top,
            where_gene="gene_name = '{gene_name}' and".format(gene_name=gene_name) if gene_name else "",
            order_by=order_by_clause
        )
        return pd.read_sql_query(sql_query, self.session.bind)


class SubclonalVariantsFigures(Figures):

    @functools.lru_cache()
    def get_top_occurring_subclonal_variants_plot(self, top, gene_name, min_vaf, order_by):

        data = SubclonalVariantsQueries(session=self.queries.session).get_top_occurring_subclonal_variants(
            top=top, gene_name=gene_name, min_vaf=min_vaf, order_by=order_by)
        fig = dcc.Markdown("""**No subclonal variants for the current selection**""")
        if data is not None and data.shape[0] > 0:

            unique_subclonal_variants = data.groupby('variant_id').agg({
                'gene_name': 'first',
                'pfam_description': 'first',
                'hgvs_p': 'first',
                'annotation_highest_impact': 'first',
                'cons_hmm_sars_cov_2': 'first',
            })
            # this fills empty conservation values
            unique_subclonal_variants["cons_hmm_sars_cov_2"].fillna(0.0, inplace=True)
            unique_subclonal_variants["cons_hmm_sars_cov_2"] = unique_subclonal_variants["cons_hmm_sars_cov_2"].transform(lambda x: round(x, 3))
            unique_subclonal_variants["median_vaf"] = data[["variant_id", "vaf"]]\
                .groupby('variant_id')["vaf"].agg(lambda x: round(float(np.median(x)), 3))
            unique_subclonal_variants["iqr_vaf"] = data[["variant_id", "vaf"]]\
                .groupby('variant_id')["vaf"].agg(lambda x: round(np.percentile(x, 0.75) - np.percentile(x, 0.25),  3))
            unique_subclonal_variants["first_observation"] = data[["variant_id", "date"]] \
                .groupby('variant_id')["date"].agg(np.min)
            unique_subclonal_variants["last_observation"] = data[["variant_id", "date"]] \
                .groupby('variant_id')["date"].agg(np.max)
            unique_subclonal_variants["count_observations"] = data[["variant_id", "date"]] \
                .groupby('variant_id')["date"].agg('count')
            unique_subclonal_variants["score"] = unique_subclonal_variants[["count_observations", "cons_hmm_sars_cov_2"]]\
                .apply(lambda x: round(np.log(x[0]) * x[1], 3), axis=1)

            unique_subclonal_variants.reset_index(inplace=True)

            # set the styles of the cells
            #styles_counts = self.discrete_background_color_bins(data, columns=included_month_colums)
            #styles_total_count = self.discrete_background_color_bins(data, columns=["total"], colors="Reds")
            #styles_frequency = self._get_table_style_by_af()
            styles_striped = [{
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }]

            if order_by == "score":
                ordered_data = unique_subclonal_variants.sort_values("score", ascending=True)
            elif order_by == "count":
                ordered_data = unique_subclonal_variants.sort_values("count_observations", ascending=False)
            elif order_by == "conservation":
                ordered_data = unique_subclonal_variants.sort_values("cons_hmm_sars_cov_2", ascending=True)
            elif order_by == "vaf":
                ordered_data = unique_subclonal_variants.sort_values("median_vaf", ascending=False)
            else:
                raise CovigatorQueryException("Not supported order by of subclonal variants")

            fig = dash_table.DataTable(
                id='top-occurring-subclonal-variants-table',
                data=ordered_data.to_dict('records'),
                columns=[
                            {"name": ["Gene"], "id": "gene_name"},
                            {"name": ["Pfam Domain"], "id": "pfam_description"},
                            {"name": ["DNA mutation"], "id": "variant_id"},
                            {"name": ["Protein mutation"], "id": "hgvs_p"},
                            {"name": ["Effect"], "id": "annotation_highest_impact"},
                            {"name": ["First observation"], "id": "first_observation"},
                            {"name": ["Last observation"], "id": "last_observation"},
                            {"name": ["Count"], "id": "count_observations"},
                            {"name": ["ConsHMM"], "id": "cons_hmm_sars_cov_2"},
                            {"name": ["Median VAF"], "id": "median_vaf"},
                            {"name": ["IQR VAF"], "id": "iqr_vaf"},
                            {"name": ["Score"], "id": "score"},
                        ],
                style_data_conditional=styles_striped,
                style_cell_conditional=[
                    {
                        'if': {'column_id': c},
                        'textAlign': 'left'
                    } for c in ['gene_name', 'pfam_description', 'variant_id', 'hgvs_p', 'annotation_highest_impact']
                ],
                style_as_list_view=True,
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
                row_selectable='single',
                css=[{'selector': '.dash-cell div.dash-cell-value',
                      'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'}],
                style_table={'overflowX': 'auto'},
                style_data={'whiteSpace': 'normal', 'height': 'auto'},
                style_cell={'maxWidth': '100px'},
            )

        return [
            fig,
            dcc.Markdown("""
            **Variants only observed as intrahost variants**
                            
            Intrahost variants are detected in the ENA dataset as variant calls with a VAF lower than 80 %; 
            all variant calls with a higher VAF are considered clonal.
            This table shows those intrahost variants that have not been observed as clonal variants. 
            Assuming that this subset of intrahost variants have evolved separately within each host some of this
            variants may recur within patients due to positive selection and hence may be in risk of becoming 
            clonal variants.
            The list of intrahost variants can be prioritised by the count of observations, 
            by the ConsHMM conservation score, by the median VAF or by a joint score using the count of observations 
            and the conservation scores (ie: *ln(count observations) x conservation score*).
                            """)
        ]