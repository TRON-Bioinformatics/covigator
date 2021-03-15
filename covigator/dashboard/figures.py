import random
import pandas as pd
import numpy as np
import plotly.express as px
from sqlalchemy import and_, func
from sqlalchemy.orm import Session
import dash_bio
from covigator.model import SampleEna, JobEna, JobStatus, Variant, Gene, VariantObservation

import json
from six.moves.urllib import request as urlreq
import dash_bio as dashbio


def get_accumulated_samples_by_country(session: Session):

    # fetch data from database
    samples = pd.read_sql(session.query(SampleEna).join(JobEna).filter(JobEna.status == JobStatus.LOADED).statement, session.bind)

    fig = None
    if samples.shape[0] > 0:
        data = _prepare_accumulated_samples(samples)
        fig = px.area(data, x="date", y="cumsum", color="country",
                      category_orders={
                          "country": list(data.sort_values("cumsum", ascending=False).country.unique())[::-1]},
                      labels={"cumsum": "num. samples", "count": "increment"},
                      title="Accumulated samples per country",
                      hover_data=["count"],
                      color_discrete_sequence=random.shuffle(px.colors.qualitative.Dark24))
        fig.update_layout(
            legend={'traceorder': 'reversed'},
            xaxis={'title': None},
            yaxis={'dtick': 2000}
        )
    return fig


def _prepare_accumulated_samples(samples: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame with columns: data, country, cumsum, count
    """

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


def get_variants_plot(session: Session, gene_name="S"):

    # reads gene annotations
    gene = session.query(Gene).filter(Gene.name == gene_name).first()
    start = int(gene.data.get("start"))
    protein_features = gene.data.get("transcripts", [])[0].get("translations", [])[0].get("protein_features")
    domains = [{"name": f.get('description'), "coord": "{}-{}".format(
        start + int(f["start"]), start + int(f["end"]))} for f in protein_features if f.get("dbname") == "Pfam"]

    # reads variants
    variants = pd.read_sql(
        session.query(VariantObservation.position, Variant.annotation, func.count(VariantObservation.position))
            .join(Variant)
            # filter out synonymous variants
            .filter(and_(Variant.gene_name == gene_name, Variant.annotation != "synonymous_variant"))
            .group_by(VariantObservation.position, Variant.annotation)
            .statement,
        session.bind)

    fig = None
    if variants.shape[0] > 0:
        # reads total number of samples and calculates frequencies
        count_samples = session.query(JobEna).filter(JobEna.status == JobStatus.LOADED).count()
        variants["af"] = variants.count_1 / count_samples
        variants["log_af"] = variants.af.transform(lambda x: np.log(x + 1))
        variants["log_count"] = variants.count_1.transform(lambda x: np.log(x))

        # TODO: do something in the data ingestion about multiple annotations on the same variant
        variants.annotation = variants.annotation.transform(lambda a: a.split("&")[0])

        mdata = {
            "x": list(variants.position.transform(lambda x: str(x))),
            "y": list(variants.log_count.transform(lambda x: round(x, 3))),
            "mutationGroups": list(variants.annotation),
            "domains": domains
        }

        fig = dash_bio.NeedlePlot(
            id='my-dashbio-needleplot',
            mutationData=mdata,
            rangeSlider=True,
            ylabel="log(count observations)",
            domainStyle={
                'displayMinorDomains': False
            }
        )
    return fig


def get_circos_plot():
    data = urlreq.urlopen(
        'https://raw.githubusercontent.com/plotly/dash-bio-docs-files/master/circos_graph_data.json').read()
    circos_graph_data = json.loads(data)

    layout = [{"id":"MN908947.3", "label": "MN908947.3", "color":"#996600", "len":29903}]

    return dashbio.Circos(
        layout=layout,
        tracks=[{
            'type': 'CHORDS',
            'data': circos_graph_data['chords']
        }],
        config={
            #'innerRadius': 40,
            #'outerRadius': 200
        }
    )
