import requests
from covigator.model import Database, EnaRun
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

    def __init__(self, tax_id: str, host_tax_id: str, database: Database):
        self.tax_id = tax_id
        assert self.tax_id is not None and self.tax_id.strip() != "", "Empty tax id"
        logger.info("Tax id {}".format(self.tax_id))
        self.host_tax_id = host_tax_id
        assert self.host_tax_id is not None and self.host_tax_id.strip() != "", "Empty host tax id"
        logger.info("Host tax id {}".format(self.host_tax_id))
        self.database = database
        assert self.database is not None, "Empty database"

        self.excluded_samples_by_host_tax_id = {}
        self.excluded_samples_by_fastq_ftp = 0
        self.excluded_samples_by_instrument_platform = {}
        self.excluded_existing = 0
        self.included = 0

    def access(self):
        offset = 0
        finished = False
        while not finished:
            list_runs = self._get_ena_runs_page(offset)
            self._process_runs(list_runs)
            if len(list_runs) < self.PAGE_SIZE:
                finished = True
            offset += len(list_runs)
        self._log_results()

    def _get_ena_runs_page(self, offset):
        return requests.get(
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

    def _process_runs(self, list_runs):

        session = self.database.get_database_session()
        included_runs = []
        try:
            for run in list_runs:
                ena_run = self._parse_ena_run(run)
                if not self._complies_with_inclusion_criteria(ena_run):
                    continue    # skips runs not complying with inclusion criteria
                if session.query(EnaRun).filter_by(run_accession=ena_run.run_accession).count() > 0:
                    self.excluded_existing += 1
                    continue    # skips runs already registered in the database
                self.included += 1
                included_runs.append(ena_run)
            if len(included_runs) > 0:
                session.add_all(included_runs)
                session.commit()
                logger.info("Added {} new runs".format(len(included_runs)))
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()

    def _parse_ena_run(self, run):
        ena_run = EnaRun(**run)
        try:
            ena_run.lat = float(ena_run.lat)
        except ValueError or TypeError:
            ena_run.lat = None
        try:
            ena_run.lon = float(ena_run.lon)
        except ValueError or TypeError:
            ena_run.lon = None
        return ena_run

    def _complies_with_inclusion_criteria(self, ena_run: EnaRun):
        included = True
        if ena_run.host_tax_id is None or ena_run.host_tax_id.strip() == "" or ena_run.host_tax_id != self.host_tax_id:
            included = False    # skips runs where the host is empty or does not match
            self.excluded_samples_by_host_tax_id[str(ena_run.host_tax_id)] = \
                self.excluded_samples_by_host_tax_id.get(str(ena_run.host_tax_id), 0) + 1
        if ena_run.fastq_ftp is None or ena_run.fastq_ftp == "":
            included = False    # skips runs without FTP URL
            self.excluded_samples_by_fastq_ftp += 1
        if ena_run.instrument_platform.upper() != "ILLUMINA":
            included = False    # skips non Illumina data
            self.excluded_samples_by_instrument_platform[str(ena_run.instrument_platform)] = \
                self.excluded_samples_by_instrument_platform.get(str(ena_run.instrument_platform), 0) + 1
        return included

    def _log_results(self):
        logger.info("Excluded existing runs = {}".format(self.excluded_existing))
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded due to empty FASTQ FTP URL runs = {}".format(self.excluded_samples_by_fastq_ftp))
        logger.info("Excluded by platform runs = {}".format(self.excluded_samples_by_instrument_platform))
        logger.info("Excluded by host if runs = {}".format(self.excluded_samples_by_host_tax_id))
