import pandas as pd
import numpy as np
from dash import dcc
from dash import html
from dash import dash_table
from logzero import logger
from sqlalchemy.orm import Session
from covigator.dashboard.figures.figures import Figures, MARGIN, TEMPLATE, PLOTLY_CONFIG, STYLES_STRIPPED, STYLE_HEADER, \
    STYLE_CELL
from covigator.database.model import SubclonalVariantObservation, VariantObservation, SampleEna
from covigator.exceptions import CovigatorQueryException
import plotly.express as px


class SubclonalVariantsQueries:

    def __init__(self, session: Session):
        self.session = session

    def get_top_occurring_subclonal_variants(self, top, gene_name, domain, min_vaf, order_by):

        if order_by == "score":
            order_by_clause = "order by score asc"
        elif order_by == "count":
            order_by_clause = "order by count_observations desc"
        elif order_by == "conservation":
            order_by_clause = "order by conservation asc"
        elif order_by == "vaf":
            order_by_clause = "order by median_vaf desc"
        else:
            raise CovigatorQueryException("Not supported order by of intrahost variants")

        # counts variants over those bins
        if domain is not None:
            where_gene = "pfam_name = '{domain}'".format(domain=domain)
        elif gene_name is not None:
            where_gene = "gene_name = '{gene_name}'".format(gene_name=gene_name)
        else:
            where_gene = "gene_name is not null"

        sql_query = """
        select variant_id, date, hgvs_p, gene_name, annotation_highest_impact, 
            cons_hmm_sars_cov_2, cons_hmm_sarbecovirus, cons_hmm_vertebrate_cov, pfam_name, pfam_description, vaf
            from {subclonal_variants_table_name} 
            where 
                vaf >= {min_vaf} and 
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
                            annotation_highest_impact != 'synonymous_variant' and
                            vaf >= {min_vaf} and 
                            {where_gene} and
                            variant_id not in (
                                select distinct variant_id 
                                from {variants_table_name} 
                                where annotation_highest_impact != 'synonymous_variant' and gene_name is not null 
                            )
                        group by variant_id
                        {order_by}
                        limit {top}) as t);
            """.format(
            min_vaf=min_vaf,
            subclonal_variants_table_name=SubclonalVariantObservation.__tablename__,
            variants_table_name=VariantObservation.__tablename__,
            top=top,
            where_gene=where_gene,
            order_by=order_by_clause
        )
        return pd.read_sql_query(sql_query, self.session.bind)

    def get_count_samples_by_library_stratgey(self, variant_id, min_vaf):
        sql_query = """
        select count(*) as count_samples, library_strategy 
        from {samples_ena_table}
        where run_accession in (
            select sample 
            from {subclonal_variant_observations_table} 
            where variant_id='{variant_id}' and vaf >= {min_vaf}
            )
        group by library_strategy
        """.format(
            variant_id=variant_id,
            min_vaf=min_vaf,
            samples_ena_table=SampleEna.__tablename__,
            subclonal_variant_observations_table=SubclonalVariantObservation.__tablename__
        )
        return pd.read_sql_query(sql_query, self.session.bind)

    def get_count_samples_by_country(self, variant_id, min_vaf):
        sql_query = """
        select count(*) as count_samples, country, collection_date as date
        from {samples_ena_table}
        where run_accession in (
            select sample 
            from {subclonal_variant_observations_table} 
            where variant_id='{variant_id}' and vaf >= {min_vaf}
            )
        group by country, collection_date
        """.format(
            variant_id=variant_id,
            min_vaf=min_vaf,
            samples_ena_table=SampleEna.__tablename__,
            subclonal_variant_observations_table=SubclonalVariantObservation.__tablename__
        )
        data = pd.read_sql_query(sql_query, self.session.bind).astype(
                {'date': 'datetime64', 'count_samples': 'float64'})

        data['cumsum'] = data.sort_values("date").groupby(['country'])['count_samples'].cumsum()
        dates = data.date[~data.date.isna()].sort_values().unique()
        countries = data[(~data.country.isna())].sort_values("cumsum", ascending=False).country.unique()

        empty_table = pd.DataFrame(
            index=pd.MultiIndex.from_product([dates, countries], names=["date", "country"]))
        empty_table["count_samples"] = 0

        # adds values into empty table
        filled_table = pd.merge(
            left=empty_table.reset_index(),
            right=data.loc[:, ["date", "country", "count_samples"]],
            on=["date", "country"],
            how='left',
            suffixes=("_x", "_y")).fillna(0)
        filled_table["count_samples"] = filled_table.count_samples_x + filled_table.count_samples_y
        filled_table['cumsum'] = filled_table.groupby(['country'])['count_samples'].cumsum()

        return filled_table

    def get_top_cooccurring_clonal_variants(self, variant_id, min_vaf, gene_name, domain):
        if domain is not None:
            where_gene = "pfam_name = '{domain}'".format(domain=domain)
        elif gene_name is not None:
            where_gene = "gene_name = '{gene_name}'".format(gene_name=gene_name)
        else:
            where_gene = "gene_name is not null"

        sql_query = """
        select count(*) as count_samples, variant_id, hgvs_p, gene_name, annotation_highest_impact
        from {variant_observations_table}
        where sample in (
            select sample 
            from {subclonal_variant_observations_table} 
            where variant_id='{variant_id}' and vaf >= {min_vaf}
            )
        and {where_gene}
        group by variant_id, hgvs_p, gene_name, annotation_highest_impact
        order by count_samples desc
        limit 10
        """.format(
            variant_id=variant_id,
            min_vaf=min_vaf,
            variant_observations_table=VariantObservation.__tablename__,
            subclonal_variant_observations_table=SubclonalVariantObservation.__tablename__,
            where_gene=where_gene
        )
        return pd.read_sql_query(sql_query, self.session.bind)


class IntrahostMutationsFigures(Figures):

    def get_top_occurring_subclonal_variants_plot(self, top, gene_name, domain, min_vaf, order_by):

        logger.debug("Getting data on top occurring intrahost mutations...")
        data = SubclonalVariantsQueries(session=self.queries.session).get_top_occurring_subclonal_variants(
            top=top, gene_name=gene_name, domain=domain, min_vaf=min_vaf, order_by=order_by)
        fig = dcc.Markdown("""**No intrahost variants for the current selection**""")
        if data is not None and data.shape[0] > 0:
            logger.debug("Preparing plot on top occurring intrahost variants...")
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

            if order_by == "score":
                ordered_data = unique_subclonal_variants.sort_values("score", ascending=True)
            elif order_by == "count":
                ordered_data = unique_subclonal_variants.sort_values("count_observations", ascending=False)
            elif order_by == "conservation":
                ordered_data = unique_subclonal_variants.sort_values("cons_hmm_sars_cov_2", ascending=True)
            elif order_by == "vaf":
                ordered_data = unique_subclonal_variants.sort_values("median_vaf", ascending=False)
            else:
                raise CovigatorQueryException("Not supported order by of intrahost variants")

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
                style_data_conditional=STYLES_STRIPPED,
                style_as_list_view=True,
                style_header=STYLE_HEADER,
                row_selectable='single',
                css=[{'selector': '.dash-cell div.dash-cell-value',
                      'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'}],
                style_table={'overflowX': 'auto'},
                style_data={'whiteSpace': 'normal', 'height': 'auto'},
                style_cell={'maxWidth': '100px'},
            )

        return [
            fig,
            html.Br(),
            html.Div(children=[
                html.Button("Download CSV", id="btn_csv3"),
                dcc.Download(id="download-dataframe-csv3"),
                dcc.Store(id="memory3", data=ordered_data.to_dict('records'))]),
            html.Br(),
            dcc.Markdown("""
            **Mutations only observed as intrahost**
                            
            The list of intrahost mutations can be prioritised by the count of observations, 
            by the ConsHMM conservation score, by the median VAF or by a joint score using the count of observations 
            and the conservation scores (ie: *ln(count observations) x conservation score*).
            The variant calls can be filtered by VAF.
            
            Select any of the mutations to explore further details.
            """)
        ]

    def get_hist_library_strategy(self, variant_id, min_vaf):

        data = SubclonalVariantsQueries(session=self.queries.session).get_count_samples_by_library_stratgey(
            variant_id=variant_id, min_vaf=min_vaf)
        fig = px.bar(data_frame=data, x="library_strategy", y="count_samples", color_discrete_sequence=["#969696"])
        fig.update_layout(
            margin=MARGIN,
            template=TEMPLATE,
            yaxis={'title': "Num. of samples"},
            xaxis={'title': None},
            legend={'title': None}
        )
        return [
            dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
            dcc.Markdown("""**Library strategies distribution for {}**""".format(variant_id))
        ]

    def get_hist_countries(self, variant_id, min_vaf):

        data = SubclonalVariantsQueries(session=self.queries.session).get_count_samples_by_country(
            variant_id=variant_id, min_vaf=min_vaf)

        countries = list(data.sort_values("cumsum", ascending=False).country.unique())
        fig = px.area(data, x="date", y="cumsum", color="country",
                      category_orders={
                          "country": countries[::-1]},
                      labels={"cumsum": "num. samples", "count_samples": "increment"},
                      hover_data=["count_samples"],
                      color_discrete_sequence=px.colors.qualitative.Light24)

        fig.update_layout(
            margin=MARGIN,
            template=TEMPLATE,
            yaxis={'title': "Num. of samples"},
            xaxis={'title': None},
            legend={'title': None}
        )
        return [
            dcc.Graph(figure=fig, config=PLOTLY_CONFIG),
            dcc.Markdown("""**Countries distribution for {}**""".format(variant_id))
        ]

    def get_cooccurring_clonal_variants(self, variant_id, min_vaf, gene_name, domain):
        data = SubclonalVariantsQueries(session=self.queries.session).get_top_cooccurring_clonal_variants(
            variant_id=variant_id, min_vaf=min_vaf, gene_name=gene_name, domain=domain)

        fig = dash_table.DataTable(
            id='top-cooccurring-clonal-variants-table',
            data=data.to_dict('records'),
            columns=[
                {"name": ["Gene"], "id": "gene_name"},
                {"name": ["DNA mutation"], "id": "variant_id"},
                {"name": ["Protein mutation"], "id": "hgvs_p"},
                {"name": ["Effect"], "id": "annotation_highest_impact"},
                {"name": ["Count samples"], "id": "count_samples"},
            ],
            style_data_conditional=STYLES_STRIPPED,
            style_as_list_view=True,
            style_header=STYLE_HEADER,
            css=[{'selector': '.dash-cell div.dash-cell-value',
                  'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'}],
            style_table={'overflowX': 'auto'},
            style_data={'whiteSpace': 'normal', 'height': 'auto'},
            style_cell=STYLE_CELL,
        )

        return [
            fig,
            dcc.Markdown("""**Top 10 co-occurring clonal mutations with the intrahost mutation {variant_id}**
                    """.format(variant_id=variant_id))
        ]
