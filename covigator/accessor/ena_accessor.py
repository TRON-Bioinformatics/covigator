from datetime import date

import pycountry
import pycountry_convert
import requests
from covigator.misc import backoff_retrier
from covigator.model import EnaRun, Job
from covigator.database import Database
from logzero import logger

NUMBER_RETRIES = 5


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
        "library_name",
        "nominal_length",
        "library_layout",
        "library_strategy",
        "library_source",
        "library_selection",
        "read_count",
        "base_count",
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
    INCLUDED_LIBRARY_STRATEGIES = [
        "WGA",
        "WGS",
        "Targeted-Capture",
        "RNA-Seq"
    ]

    def __init__(self, tax_id: str, host_tax_id: str, database: Database, maximum=None):
        self.tax_id = tax_id
        assert self.tax_id is not None and self.tax_id.strip() != "", "Empty tax id"
        logger.info("Tax id {}".format(self.tax_id))
        self.host_tax_id = host_tax_id
        assert self.host_tax_id is not None and self.host_tax_id.strip() != "", "Empty host tax id"
        logger.info("Host tax id {}".format(self.host_tax_id))
        self.database = database
        assert self.database is not None, "Empty database"
        self.maximum = maximum

        self.excluded_samples_by_host_tax_id = {}
        self.excluded_samples_by_fastq_ftp = 0
        self.excluded_samples_by_instrument_platform = {}
        self.excluded_samples_by_library_strategy = {}
        self.excluded_existing = 0
        self.included = 0
        self.excluded = 0

        # this ensures there is a retry mechanism in place with a limited number of retries
        self.get_with_retries = backoff_retrier.wrapper(requests.get, NUMBER_RETRIES)

    def access(self):
        offset = 0
        finished = False
        while not finished:
            list_runs = self._get_ena_runs_page(offset)
            self._process_runs(list_runs)
            # finishes when no more data or when test parameter maximum is reached
            if len(list_runs) < self.PAGE_SIZE or (self.maximum is not None and self.included >= self.maximum):
                finished = True
            offset += len(list_runs)
        self._log_results()

    def _get_ena_runs_page(self, offset):
        return self.get_with_retries(
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
        included_jobs = []
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
                included_jobs.append(Job(run_accession=ena_run.run_accession))
            if len(included_runs) > 0:
                session.add_all(included_runs)
                session.commit()    # we need two commits to maintain integrity of FKs
                session.add_all(included_jobs)
                session.commit()
                logger.info("Added {} new runs".format(len(included_runs)))
        except Exception as e:
            logger.exception(e)
            session.rollback()
        finally:
            session.close()

    def _parse_country(self, ena_run):
        ena_run.country_raw = ena_run.country
        if ena_run.country_raw is not None and ena_run.country_raw.strip() == "":
            ena_run.country_raw = None
        try:
            match = pycountry.countries.search_fuzzy(ena_run.country_raw.split(":")[0])[0]
            ena_run.country = match.name
            ena_run.country_alpha_2 = match.alpha_2
            ena_run.country_alpha_3 = match.alpha_3
            ena_run.continent_alpha_2 = pycountry_convert.country_alpha2_to_continent_code(ena_run.country_alpha_2)
            ena_run.continent = pycountry_convert.convert_continent_code_to_continent_name(ena_run.continent_alpha_2)
        except (LookupError, AttributeError):
            #logger.error("Error parsing country {}. Filling with NAs".format(ena_run.country))
            ena_run.country = "Not available"
            ena_run.country_alpha_2 = "None"    # don't use NA as that is the id of Namibia
            ena_run.country_alpha_3 = "None"
            ena_run.continent_alpha_2 = "None"
            ena_run.continent = "None"

    def _parse_dates(self, ena_run):
        ena_run.collection_date = self._parse_abstract(ena_run.collection_date, date.fromisoformat)
        ena_run.first_created = self._parse_abstract(ena_run.first_created, date.fromisoformat)

    def _parse_ena_run(self, run):
        ena_run = EnaRun(**run)
        self._parse_country(ena_run)
        self._parse_dates(ena_run)
        self._parse_numeric_fields(ena_run)
        fastqs = ena_run.get_fastqs_ftp()
        # annotates with the number of FASTQ files, this is useful as we hold the FASTQs in a single string
        ena_run.num_fastqs = 0 if fastqs is None or fastqs == [""] else len(fastqs)
        return ena_run

    def _parse_numeric_fields(self, ena_run):
        ena_run.nominal_length = self._parse_abstract(ena_run.nominal_length, int)
        ena_run.read_count = self._parse_abstract(ena_run.read_count, int)
        ena_run.base_count = self._parse_abstract(ena_run.base_count, int)
        ena_run.lat = self._parse_abstract(ena_run.lat, float)
        ena_run.lon = self._parse_abstract(ena_run.lon, float)

    def _parse_abstract(self, value, type):
        try:
            value = type(value)
        except (ValueError, TypeError):
            value = None
        return value

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
        if ena_run.library_strategy not in self.INCLUDED_LIBRARY_STRATEGIES:
            included = False  # skips not included library strategies
            self.excluded_samples_by_library_strategy[str(ena_run.library_strategy)] = \
                self.excluded_samples_by_library_strategy.get(str(ena_run.library_strategy), 0) + 1
        if not included:
            self.excluded += 1
        return included

    def _log_results(self):
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded already existing runs = {}".format(self.excluded_existing))
        logger.info("Total excluded runs by selection criteria = {}".format(self.excluded))
        logger.info("Excluded due to empty FASTQ FTP URL runs = {}".format(self.excluded_samples_by_fastq_ftp))
        logger.info("Excluded by platform runs = {}".format(self.excluded_samples_by_instrument_platform))
        logger.info("Excluded by host if runs = {}".format(self.excluded_samples_by_host_tax_id))
        logger.info("Excluded by library strategy = {}".format(self.excluded_samples_by_library_strategy))
