import csv
from datetime import date, datetime
import json

import pycountry
import pycountry_convert
import requests
from sqlalchemy.orm import Session

from COVIGATOR import ENV_META_GISAID

from covigator.misc import backoff_retrier
from covigator.database.model import SampleGisaid, JobGisaid, Sample, DataSource, Log, CovigatorModule
from covigator.database.database import Database
from logzero import logger


class GisaidAccessor:
    GISAID_METADATA_FILE = os.getenv(ENV_COVIGATOR_REF_GISAID, "/scratch/info/projects/SARS-CoV-2/gisaid/gisaid_hcov-19_2020_10_02_11_ST_corrected_v2.tsv")
    
    GISAID_FIELDS = [
        # data on run
        "virus_name",
        "collection_date",
        "instrument_platform",
        "instrument_model",
        "assembly_method",
        # data on host
        "host",
        "host_body_site",
        # geographical data
        "lat",
        "lon",
        "country"
    ]

    INCLUDED_INSTRUMENT_PLATFORMS = [
        "ILLUMINA",
        "ION TORRENT",
        "SANGER",
        "MGISEQ",
        "NANOPORE",
        "BGISEQ",
        "BIOELECTRONSEQ",
        "PACBIO",
        "THERMOFISHER",
        "UNKNOWN"
    ]

    # Not normalized within GISAID data; too many; will not be used for now
    INCLUDED_ASSEMBLY_METHODS = [
        "ARTICv3",
        "BBMap",
        "CLC Genomics workbench",
        "DNAStar",
        "Seattle Flu Assembly Pipeline",
        "Geneious Prime"
    ]

    def __init__(self, tax_id: str, host_tax_id: str, database: Database, maximum=None):
        self.start_time = datetime.now()
        self.has_error = False
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
        self.excluded_samples_by_instrument_platform = {}
        self.excluded_samples_by_assembly_method = {}
        self.excluded_existing = 0
        self.included = 0
        self.excluded = 0

        # this ensures there is a retry mechanism in place with a limited number of retries
        #self.get_with_retries = backoff_retrier.wrapper(requests.get, NUMBER_RETRIES)

    def access(self):
        session = self.database.get_database_session()
        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # it assumes there will be no repetitions
        existing_sample_ids = [value for value, in session.query(SampleGisaid.run_accession).all()]
        try:
            list_runs = self._get_gisaid_runs()
            #print(list_runs)
            logger.info("Read file of {} GISAID samples".format(len(list_runs)))
            self._process_runs(list_runs, existing_sample_ids, session)

        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.has_error = True
        finally:
            self._write_execution_log(session)
            session.close()
            self._log_results()

    def _get_gisaid_runs(self):
        # Exemplary columns to be read
        # 'Virus name': 'hCoV-19/USA/MA-MGH-00523/2020', 
        # 'Accession ID': 'EPI_ISL_460262', 
        # 'Collection date': '2020-03-18', 
        # 'Location': 'North America / USA / Massachusetts', 
        # 'Host': 'Human', 
        # 'Passage': 'Original', 
        # 'Specimen': 'oronasopharynx', 
        # 'Additional host information': '', 
        # 'Sequencing technology': 'Illumina NovaSeq', 
        # 'Assembly method': '', 
        # 'Comment': '', 
        # 'Comment type': '', 
        # 'Lineage': 'B.1', 
        # 'Clade': 'GH'

        results = []
        with open(self.GISAID_METADATA_FILE) as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            #{'run_accession': 'EPI_ISL_463264', 
            # 'virus_name': 'hCoV-19/Canada/BC_27280425/2020', 
            # 'collection_date': '2020-03', 
            # 'host': 'Human', 
            # 'host_body_site': 'NA', 
            # 'instrument_platform': 'nanopore', 
            # 'instrument_model': 'oxford nanopore', 
            # 'country': 'Canada'}
            for row in reader:
                fields = {}
                fields["host_tax_id"] = self.host_tax_id
                fields["run_accession"] = row["accession_id"]
                fields["virus_name"] = row["virus_name"]
                fields["collection_date"] = row["collection_date"]
                fields["host"] = row["host"]
                fields["host_body_site"] = row["specimen"]
                fields["instrument_platform"] = row["sequencing_class"]
                #if "ILLUMINA" in row["sequencing_technology"].upper():
                #    fields["instrument_platform"] = "ILLUMINA"
                #else:
                #    fields["instrument_platform"] = "unknown"
                fields["instrument_model"] = row["sequencing_technology"]
                #fields["assembly_method"] = row["assembly_method"]
                fields["country"] = row["location"].split("/")[1].strip()
                results.append(fields)
        return results

    def _process_runs(self, list_runs, existing_sample_ids, session: Session):

        included_samples = []
        included_samples_gisaid = []
        included_jobs = []
        for run in list_runs:
            if run.get("run_accession") in existing_sample_ids:
                self.excluded_existing += 1
                continue    # skips runs already registered in the database
            if not self._complies_with_inclusion_criteria(run):
                continue    # skips runs not complying with inclusion criteria
            # NOTE: this parse operation is costly
            sample_gisaid = self._parse_gisaid_run(run)
            sample = self._build_sample(sample_gisaid)
            self.included += 1
            included_samples_gisaid.append(sample_gisaid)
            included_samples.append(sample)
            included_jobs.append(JobGisaid(run_accession=sample_gisaid.run_accession))
        if len(included_samples) > 0:
            session.add_all(included_samples_gisaid)
            session.commit()
            session.add_all(included_samples)
            session.add_all(included_jobs)
            session.commit()
            logger.info("Added {} new GISAID samples".format(len(included_samples)))
        logger.info("Processed {} GISAID samples".format(len(list_runs)))

    def _build_sample(self, sample_gisaid):
        return Sample(
            id=sample_gisaid.run_accession,
            source=DataSource.GISAID,
            gisaid_id=sample_gisaid.run_accession
        )

    def _parse_country(self, gisaid_run):
        gisaid_run.country_raw = gisaid_run.country
        if gisaid_run.country_raw is not None and gisaid_run.country_raw.strip() == "":
            gisaid_run.country_raw = None
        try:
            match = pycountry.countries.search_fuzzy(gisaid_run.country_raw.split("/")[0])[0]
            gisaid_run.country = match.name
            gisaid_run.country_alpha_2 = match.alpha_2
            gisaid_run.country_alpha_3 = match.alpha_3
            gisaid_run.continent_alpha_2 = pycountry_convert.country_alpha2_to_continent_code(gisaid_run.country_alpha_2)
            gisaid_run.continent = pycountry_convert.convert_continent_code_to_continent_name(gisaid_run.continent_alpha_2)
        except (LookupError, AttributeError):
            #logger.error("Error parsing country {}. Filling with NAs".format(gisaid_run.country))
            gisaid_run.country = "Not available"
            gisaid_run.country_alpha_2 = "None"    # don't use NA as that is the id of Namibia
            gisaid_run.country_alpha_3 = "None"
            gisaid_run.continent_alpha_2 = "None"
            gisaid_run.continent = "None"

    def _parse_dates(self, gisaid_run):
        gisaid_run.collection_date = self._parse_abstract(gisaid_run.collection_date, date.fromisoformat)
        gisaid_run.first_created = self._parse_abstract(gisaid_run.first_created, date.fromisoformat)

    def _parse_gisaid_run(self, run):
        gisaid_run = SampleGisaid(**run)
        self._parse_country(gisaid_run)
        self._parse_dates(gisaid_run)
        self._parse_numeric_fields(gisaid_run)
        return gisaid_run

    def _parse_numeric_fields(self, gisaid_run):
        gisaid_run.lat = self._parse_abstract(gisaid_run.lat, float)
        gisaid_run.lon = self._parse_abstract(gisaid_run.lon, float)

    def _parse_abstract(self, value, type):
        try:
            value = type(value)
        except (ValueError, TypeError):
            value = None
        return value

    def _complies_with_inclusion_criteria(self, gisaid_run: dict):
        # NOTE: this uses the original dictionary instead of the parsed SampleGisaid class for performance reasons
        included = True
        host_tax_id = gisaid_run.get("host_tax_id")
        if host_tax_id is None or host_tax_id.strip() == "" or host_tax_id != self.host_tax_id:
            included = False    # skips runs where the host is empty or does not match
            self.excluded_samples_by_host_tax_id[str(host_tax_id)] = \
                self.excluded_samples_by_host_tax_id.get(str(host_tax_id), 0) + 1
        # too messy within GISAID metadata; requires normalization; maybe skip for now
        instrument_platform = gisaid_run.get("instrument_platform")
        if instrument_platform.upper() not in self.INCLUDED_INSTRUMENT_PLATFORMS:
            included = False    # skips non Illumina data
            self.excluded_samples_by_instrument_platform[str(instrument_platform)] = \
                self.excluded_samples_by_instrument_platform.get(str(instrument_platform), 0) + 1
        # too messy within GISAID metadata; requires normalization; maybe skip for now
        #assembly_method = gisaid_run.get("assembly_method")
        #if assembly_method not in self.INCLUDED_ASSEMBLY_METHODS:
        #    included = False  # skips not included library strategies
        #    self.excluded_samples_by_assembly_method[str(assembly_method)] = \
        #        self.excluded_samples_by_assembly_method.get(str(assembly_method), 0) + 1
        if not included:
            self.excluded += 1
        return included

    def _log_results(self):
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded already existing runs = {}".format(self.excluded_existing))
        logger.info("Total excluded runs by selection criteria = {}".format(self.excluded))
        # logger.info("Excluded due to empty FASTQ FTP URL runs = {}".format(self.excluded_samples_by_fastq_ftp))
        logger.info("Excluded by platform runs = {}".format(self.excluded_samples_by_instrument_platform))
        logger.info("Excluded by host if runs = {}".format(self.excluded_samples_by_host_tax_id))
        logger.info("Excluded by assembly method = {}".format(self.excluded_samples_by_assembly_method))

    def _write_execution_log(self, session: Session):
        end_time = datetime.now()
        session.add(Log(
            start=self.start_time,
            end=end_time,
            source=DataSource.GISAID,
            module=CovigatorModule.ACCESSOR,
            has_error=self.has_error,
            data={
                "included": self.included,
                "excluded": {
                    "existing": self.excluded_existing,
                    "excluded_by_criteria": self.excluded,
                    #"platform": self.excluded_samples_by_instrument_platform,
                    #"assembly_method": self.excluded_samples_by_assembly_method,
                    "host": self.excluded_samples_by_host_tax_id
                }
            }
        ))
        session.commit()
