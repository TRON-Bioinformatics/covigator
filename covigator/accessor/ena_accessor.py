from datetime import datetime
from json import JSONDecodeError
from requests import Response
from sqlalchemy.orm import Session
import covigator
from covigator.accessor.abstract_accessor import AbstractAccessor, _parse_abstract
from covigator.exceptions import CovigatorExcludedSampleTooEarlyDateException
from covigator.database.model import SampleEna, DataSource, Log, CovigatorModule
from covigator.database.database import Database
from logzero import logger

BATCH_SIZE = 1000


def _parse_numeric_fields(ena_run):
    ena_run.nominal_length = _parse_abstract(ena_run.nominal_length, int)
    ena_run.read_count = _parse_abstract(ena_run.read_count, int)
    ena_run.base_count = _parse_abstract(ena_run.base_count, int)
    ena_run.lat = _parse_abstract(ena_run.lat, float)
    ena_run.lon = _parse_abstract(ena_run.lon, float)


class EnaAccessor(AbstractAccessor):

    ENA_API_URL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
    # see https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=json for all possible fields
    ENA_FIELDS = [
        # data on run
        "scientific_name",
        "sample_accession",
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
        "tax_id",
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
        "Targeted-Capture"
    ]

    def __init__(self, tax_id: str, host_tax_id: str, database: Database, maximum=None,
                 disable_minimum_date: bool=False):

        super().__init__()
        logger.info("Initialising ENA accessor")
        self.start_time = datetime.now()
        self.has_error = False
        self.error_message = None
        self.tax_id = tax_id
        assert self.tax_id is not None and self.tax_id.strip() != "", "Empty tax id"
        logger.info("Tax id {}".format(self.tax_id))
        if host_tax_id is not None and host_tax_id.strip() != "":
            self.host_tax_id = host_tax_id
            self.host_tax_id_filter = True
            logger.info("Host tax id {}".format(self.host_tax_id))
        else:
            self.host_tax_id = None
            self.host_tax_id_filter = False
            logger.info("Empty host tax id. Disabling host id filter criteria")
        self.database = database
        assert self.database is not None, "Empty database"
        self.maximum = maximum
        self.disable_minimum_date = disable_minimum_date
        if self.disable_minimum_date:
            logger.info("Disabling minimum date filter criteria")

        self.excluded_samples_by_host_tax_id = {}
        self.excluded_samples_by_fastq_ftp = 0
        self.excluded_samples_by_instrument_platform = {}
        self.excluded_samples_by_library_strategy = {}
        self.excluded_existing = 0
        self.included = 0
        self.excluded = 0
        self.excluded_by_date = 0

    def access(self):
        logger.info("Starting ENA accessor")
        session = self.database.get_database_session()
        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # it assumes there will be no repetitions
        existing_sample_ids = [value for value, in session.query(SampleEna.run_accession).all()]
        try:
            logger.info("Reading...")
            list_runs = self._get_ena_runs_page()
            logger.info("Processing {} ENA samples...".format(len(list_runs)))
            self._process_runs(list_runs, existing_sample_ids, session)
            logger.info("All samples processed!")
        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.has_error = True
            self.error_message = str(e)
        finally:
            self._write_execution_log(session)
            session.close()
            self._log_results()
            logger.info("Finished ENA accessor")

    def _get_ena_runs_page(self):
        # as communicated by ENA support we use limit=0 and offset=0 to get all records in one query
        response: Response = self.get_with_retries(
            "{url_base}/search?result=read_run&"
            "query=tax_eq({tax_id})&"
            "limit=0&"
            "offset=0&"
            "fields={fields}&"
            "format=json".format(url_base=self.ENA_API_URL_BASE, tax_id=self.tax_id, fields=",".join(self.ENA_FIELDS)))
        try:
            json = response.json()
        except JSONDecodeError as e:
            logger.exception(e)
            raise e
        return json

    def _process_runs(self, list_runs, existing_sample_ids, session: Session):

        included_samples = []
        for run in list_runs:
            if isinstance(run, dict):
                if run.get("run_accession") in existing_sample_ids:
                    self.excluded_existing += 1
                    continue    # skips runs already registered in the database
                if not self._complies_with_inclusion_criteria(run):
                    continue    # skips runs not complying with inclusion criteria
                # NOTE: this parse operation is costly
                # Remove tax_id field introduced with 2023-05-02 API update
                if "tax_id" in run:
                    run.pop("tax_id")
                try:
                    sample_ena = self._parse_ena_run(run)
                    self.included += 1
                    included_samples.append(sample_ena)
                except CovigatorExcludedSampleTooEarlyDateException:
                    logger.error("Excluded sample due to too early date")
                    self.excluded_by_date += 1
                    self.excluded += 1
            else:
                logger.error("Run from ENA without the expected format")

            if len(included_samples) >= BATCH_SIZE:
                session.add_all(included_samples)
                session.commit()
                included_samples = []

        if len(included_samples) > 0:
            session.add_all(included_samples)
            session.commit()
            logger.info("Added {} new ENA samples".format(len(included_samples)))
        logger.info("Processed {} ENA samples".format(len(list_runs)))

    def _parse_ena_run(self, run):
        sample = SampleEna(**run)
        self._parse_country(sample)
        self._parse_dates(sample, self.disable_minimum_date)
        _parse_numeric_fields(sample)
        fastqs = sample.get_fastqs_ftp()
        # annotates with the number of FASTQ files, this is useful as we hold the FASTQs in a single string
        sample.num_fastqs = 0 if fastqs is None or fastqs == [""] else len(fastqs)
        sample.covigator_accessor_version = covigator.VERSION
        return sample

    def _complies_with_inclusion_criteria(self, ena_run: dict):
        # NOTE: this uses the original dictionary instead of the parsed SampleEna class for performance reasons
        included = True
        # Skip host id filter if data is not available for selected Virus
        if self.host_tax_id_filter:
            host_tax_id = ena_run.get("host_tax_id")
            if host_tax_id is None or host_tax_id.strip() == "" or host_tax_id != self.host_tax_id:
                included = False    # skips runs where the host is empty or does not match
                self.excluded_samples_by_host_tax_id[str(host_tax_id)] = \
                    self.excluded_samples_by_host_tax_id.get(str(host_tax_id), 0) + 1
        fastq_ftp = ena_run.get("fastq_ftp")
        if fastq_ftp is None or fastq_ftp == "":
            included = False    # skips runs without FTP URL
            self.excluded_samples_by_fastq_ftp += 1
        instrument_platform = ena_run.get("instrument_platform")
        if instrument_platform.upper() != "ILLUMINA":
            included = False    # skips non Illumina data
            self.excluded_samples_by_instrument_platform[str(instrument_platform)] = \
                self.excluded_samples_by_instrument_platform.get(str(instrument_platform), 0) + 1
        library_strategy = ena_run.get("library_strategy")
        if library_strategy not in self.INCLUDED_LIBRARY_STRATEGIES:
            included = False  # skips not included library strategies
            self.excluded_samples_by_library_strategy[str(library_strategy)] = \
                self.excluded_samples_by_library_strategy.get(str(library_strategy), 0) + 1
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

    def _write_execution_log(self, session: Session):
        end_time = datetime.now()
        session.add(Log(
            start=self.start_time,
            end=end_time,
            source=DataSource.ENA,
            module=CovigatorModule.ACCESSOR,
            has_error=self.has_error,
            error_message=self.error_message,
            processed=self.included,
            data={
                "included": self.included,
                "excluded": {
                    "existing": self.excluded_existing,
                    "excluded_by_criteria": self.excluded,
                    "excluded_by_date": self.excluded_by_date,
                    "missing_fastq": self.excluded_samples_by_fastq_ftp,
                    "platform": self.excluded_samples_by_instrument_platform,
                    "library_strategy": self.excluded_samples_by_library_strategy,
                    "host": self.excluded_samples_by_host_tax_id
                }
            }
        ))
        session.commit()
