import requests
from covigator.model import EnaRun
from dask.distributed import Client
from logzero import logger


class EnaDownloader:

    ENA_API_URL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
    PAGE_SIZE = 100

    def __init__(self, tax_id: str, host_tax_id: str, dask_client: Client, data_folder: str):
        self.tax_id = tax_id
        self.host_tax_id = host_tax_id
        self.dask_client = dask_client
        self.data_folder = data_folder

    def download(self):
        offset = 0
        finished = False
        while not finished:
            list_runs = requests.get(
                "{url_base}/search?result=read_run&"
                "query=tax_eq({tax_id})&"
                "limit={page_size}&"
                "offset={offset}&"
                "fields=fastq_ftp,host_tax_id,host_sex,lat,lon,country&"    # TODO: review which fields we want. All?
                "format=json".format(
                    url_base=self.ENA_API_URL_BASE, tax_id=self.tax_id, page_size=self.PAGE_SIZE, offset=offset)).json()
            self._process_runs(list_runs)
            if len(list_runs) < self.PAGE_SIZE:
                finished = True
            offset += len(list_runs)

    def _process_runs(self, list_runs):
        for run in list_runs:
            ena_run = EnaRun(**run)
            if ena_run.host_tax_id != self.host_tax_id:
                continue    # skips runs where the host does not match
            if ena_run.fastq_ftp is None or ena_run.fastq_ftp == "":
                continue  # skips runs without FTP URL
            # TODO: check if the object exists already in the DB in which case skip
            # TODO: are there any other exclusion criteria? for instance sequencing technology
            self.dask_client.submit(EnaDownloader._download_run, ena_run)

    @staticmethod
    def _download_run(ena_run: EnaRun):
        # TODO: download the file somewhere
        logger.info(ena_run.run_accession)
        logger.info(ena_run.fastq_ftp)
        # TODO: update the entry in the database

