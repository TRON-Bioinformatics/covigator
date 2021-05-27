import base64
from datetime import date, datetime
import pycountry
import pycountry_convert
from sqlalchemy.orm import Session
from typing import List
from covigator.database.model import SampleGisaid, JobGisaid, Sample, DataSource, Log, CovigatorModule
from covigator.database.database import Database
from logzero import logger
from Bio import SeqIO
import zlib


class GisaidAccessor:

    def __init__(self, input_fasta: str, input_metadata: str, database: Database):
        self.start_time = datetime.now()
        self.input_fasta = input_fasta
        self.input_metadata = input_metadata
        self.has_error = False

        self.database = database
        assert self.database is not None, "Empty database"

        self.excluded_by_host = 0
        self.excluded_existing = 0
        self.included = 0

    def access(self):
        session = self.database.get_database_session()
        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # it assumes there will be no repetitions
        existing_sample_ids = [value for value, in session.query(SampleGisaid.run_accession).all()]
        try:
            samples = self._get_gisaid_samples(existing_sample_ids)
            logger.info("Read file of {} GISAID samples".format(len(samples)))
            self._process_runs(samples, session)
        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.has_error = True
        finally:
            self._write_execution_log(session)
            session.close()
            self._log_results()

    def _get_gisaid_samples(self, existing_sample_ids) -> List[SampleGisaid]:
        meta_dict = {}
        with open(self.input_metadata) as metadata:
            for line in metadata:
                elements = line.rstrip().split("\t")
                virus_name = elements[0]
                virus_type = elements[1]
                accession_id = elements[2]
                collection_date = elements[3]
                location = elements[4]
                host = elements[7]
                meta_dict[accession_id] = (virus_type, accession_id, collection_date, location, host)


        results = {}
        for record in SeqIO.parse(self.input_fasta, "fasta"):
            fields = record.description.split("|")
            identifier = fields[3]
            if not identifier in meta_dict:
                continue
            metadata = meta_dict[identifier]
            location_list = metadata[3].split("/")
            if identifier in existing_sample_ids:
                # we skip samples already in the database
                self.excluded_existing += 1
                continue
            host = metadata[4].lower().strip()
            if host != "human":
                self.excluded_by_host += 1
                continue
            if identifier in results:
                sample_gisaid = results.get(identifier)
            else:
                sample_gisaid = SampleGisaid(
                    run_accession=identifier,
                    date=metadata[2],
                    host_tax_id=None,
                    host=host,
                    country_raw=location_list[1],
                    region=location_list[2],
                    country=None,
                    country_alpha_2=None,
                    country_alpha_3=None,
                    continent=location_list[0],
                    continent_alpha_2=None,
                    site=None,
                    site2=None,
                    sequence={}
                )

                self._parse_country(sample_gisaid)
                self._parse_dates(sample_gisaid)
                results[sample_gisaid.run_accession] = sample_gisaid
            # stores the sequence in dictionary with the gene name
            sample_gisaid.sequence[fields[0]] = self.compress_sequence(record.seq)
        return list(results.values())

    @staticmethod
    def compress_sequence(sequence):
        return base64.b64encode(zlib.compress(sequence.encode('utf-8'))).decode()

    def _process_runs(self, list_samples: List[SampleGisaid], session: Session):

        included_samples = []
        included_samples_gisaid = []
        included_jobs = []
        for sample_gisaid in list_samples:
            # NOTE: this parse operation is costly
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
        logger.info("Processed {} GISAID samples".format(len(list_samples)))

    def _build_sample(self, sample_gisaid):
        return Sample(
            id=sample_gisaid.run_accession,
            source=DataSource.GISAID,
            gisaid_id=sample_gisaid.run_accession
        )

    def _parse_country(self, gisaid_sample: SampleGisaid):
        if gisaid_sample.country_raw is not None and gisaid_sample.country_raw.strip() == "":
            gisaid_sample.country_raw = None
        try:
            match = pycountry.countries.search_fuzzy(gisaid_sample.country_raw)[0]
            gisaid_sample.country = match.name
            gisaid_sample.country_alpha_2 = match.alpha_2
            gisaid_sample.country_alpha_3 = match.alpha_3
            gisaid_sample.continent_alpha_2 = pycountry_convert.country_alpha2_to_continent_code(gisaid_sample.country_alpha_2)
            gisaid_sample.continent = pycountry_convert.convert_continent_code_to_continent_name(gisaid_sample.continent_alpha_2)
        except (LookupError, AttributeError):
            #logger.error("Error parsing country {}. Filling with NAs".format(gisaid_run.country))
            gisaid_sample.country = "Not available"
            gisaid_sample.country_alpha_2 = "None"    # don't use NA as that is the id of Namibia
            gisaid_sample.country_alpha_3 = "None"
            gisaid_sample.continent_alpha_2 = "None"
            gisaid_sample.continent = "None"

    def _parse_dates(self, gisaid_run: SampleGisaid):
        gisaid_run.date = self._parse_abstract(gisaid_run.date, date.fromisoformat)

    def _parse_abstract(self, value, type):
        try:
            value = type(value)
        except (ValueError, TypeError):
            value = None
        return value

    def _log_results(self):
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded already existing runs = {}".format(self.excluded_existing))
        logger.info("Excluded non human host = {}".format(self.excluded_by_host))

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
                    "excluded_by_host": self.excluded_by_host
                }
            }
        ))
        session.commit()
