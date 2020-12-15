import requests
from covigator.model import EnaRun
from logzero import logger


class EnaAccessor:

    ENA_API_URL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
    PAGE_SIZE = 1000
    # see https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=json for all possible fields
    ENA_FIELDS = [
        # data on run
        "scientific_name",
        "study_accession",
        "experiment_accession",
        "first_created",
        "collection_date",
        "instrument_platform",
        "instrument_model",
        "sample_collection",
        "sequencing_method",
        "center_name",
        # FASTQs
        "fastq_ftp",
        "fastq_md5",
        # data on host
        "host_tax_id",
        "host_sex",
        "host_body_site",
        "host_gravidity",
        "host_phenotype",
        "host_genotype",
        # geographical data
        "lat",
        "lon",
        "country"
    ]

    def __init__(self, tax_id: str, host_tax_id: str):
        self.tax_id = tax_id
        assert self.tax_id is not None and self.tax_id.strip() != "", "Empty tax id"
        self.host_tax_id = host_tax_id
        assert self.host_tax_id is not None and self.host_tax_id.strip() != "", "Empty host tax id"

    def access(self):
        offset = 0
        finished = False
        while not finished:
            list_runs = requests.get(
                "{url_base}/search?result=read_run&"
                "query=tax_eq({tax_id})&"
                "limit={page_size}&"
                "offset={offset}&"
                "fields={fields}&"
                "format=json".format(
                    url_base=self.ENA_API_URL_BASE,
                    tax_id=self.tax_id,
                    page_size=self.PAGE_SIZE,
                    offset=offset,
                    fields=",".join(self.ENA_FIELDS)
                )).json()
            self._process_runs(list_runs)
            if len(list_runs) < self.PAGE_SIZE:
                finished = True
            offset += len(list_runs)

    def _process_runs(self, list_runs):
        for run in list_runs:
            ena_run = EnaRun(**run)
            if not self._complies_with_inclusion_criteria(ena_run):
                continue
            self._store_run(ena_run)

    def _complies_with_inclusion_criteria(self, ena_run: EnaRun):
        included = True
        if ena_run.host_tax_id is None or ena_run.host_tax_id.strip() == "" or ena_run.host_tax_id != self.host_tax_id:
            included = False    # skips runs where the host is empty or does not match
        if ena_run.fastq_ftp is None or ena_run.fastq_ftp == "":
            included = False    # skips runs without FTP URL
        if ena_run.instrument_platform.upper() != "ILLUMINA":
            included = False    # skips non Illumina data
        # TODO: are there any other exclusion criteria? for instance sequencing technology
        return included

    def _store_run(self, ena_run: EnaRun):
        # TODO: check if the object exists already in the DB in which case skip
        # TODO: update the entry in the database
        logger.info(ena_run.__dict__)
        #logger.info(ena_run.fastq_ftp)
        #logger.info(ena_run.host_tax_id)


