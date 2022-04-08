from datetime import date, datetime
from sqlalchemy.orm import Session

import covigator
from covigator.accessor import MINIMUM_DATE
from covigator.database.model import SampleGisaid, DataSource, Log, CovigatorModule, JobStatus
from covigator.database.database import Database
from logzero import logger
from Bio import SeqIO
import time

from covigator.exceptions import CovigatorExcludedSampleTooEarlyDateException
from covigator.misc.compression import compress_sequence
from covigator.misc.country_parser import CountryParser

THRESHOLD_NON_VALID_BASES_RATIO = 0.2
THRESHOLD_GENOME_COVERAGE = 0.2
CHROMOSOME = "MN908947.3"
GENOME_LENGTH = 29903
BATCH_SIZE = 1000


class GisaidAccessor:

    def __init__(self, input_fasta: str, input_metadata: str, database: Database):
        self.start_time = datetime.now()
        self.input_fasta = input_fasta
        self.input_metadata = input_metadata
        self.has_error = False

        self.database = database
        assert self.database is not None, "Empty database"

        self.excluded_by_host = 0
        self.excluded_by_horizontal_coverage = 0
        self.excluded_by_non_valid_bases_ratio = 0
        self.excluded_existing = 0
        self.excluded_by_date = 0
        self.included = 0
        self.country_parser = CountryParser()

    def access(self):
        session = self.database.get_database_session()
        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # it assumes there will be no repetitions
        existing_sample_ids = [value for value, in session.query(SampleGisaid.run_accession).all()]
        try:
            num_samples = self.process(existing_sample_ids, session)
            logger.info("Read file of {} GISAID samples".format(num_samples))
        except Exception as e:
            logger.exception(e)
            session.rollback()
            self.has_error = True
        finally:
            self._write_execution_log(session)
            session.close()
            self._log_results()
            logger.info("Finished GISAID accessor")

    def process(self, existing_sample_ids, session) -> int:
        meta_dict = {}
        logger.info("Reading metadata file...")
        with open(self.input_metadata) as metadata:
            for line in metadata:
                elements = line.rstrip().split("\t")
                virus_name = elements[0]
                virus_type = elements[1]
                accession_id = elements[2]
                collection_date = elements[3]
                location = elements[4]
                host = elements[7]
                meta_dict[virus_name] = (virus_type, accession_id, collection_date, location, host)
        logger.info("Metadata (num_samples={}) loaded.".format(len(meta_dict)))

        num_samples = 0
        gisaid_samples = set()
        total_time = 0
        logger.info("Reading FASTA...")
        samples_gisaid = []
        for record in SeqIO.parse(self.input_fasta, "fasta"):
            start = time.time()

            fields = record.description.split("|")
            identifier = fields[0]

            if identifier in gisaid_samples:
                continue
            gisaid_samples.add(identifier)

            if identifier in existing_sample_ids:
                num_samples += 1
                # we skip samples already in the database
                self.excluded_existing += 1
                continue

            metadata = meta_dict.get(identifier)
            if not metadata:
                continue
            
            location_list = metadata[3].split("/")
            country_raw = "None"
            try:
                country_raw = location_list[1].strip()
            except:
                pass
            region = "None"
            try:
                region = location_list[2].strip()
            except:
                pass

            host = metadata[4].lower().strip()
            if host != "human":
                num_samples += 1
                self.excluded_by_host += 1
                continue

            sequence_length = len(record.seq)
            if float(sequence_length) / GENOME_LENGTH < THRESHOLD_GENOME_COVERAGE:
                num_samples += 1
                self.excluded_by_horizontal_coverage += 1
                continue

            count_n_bases = record.seq.count("N")
            count_ambiguous_bases = sum([record.seq.count(b) for b in "RYWSMKHBVD"])
            if float(count_n_bases + count_ambiguous_bases) / sequence_length > THRESHOLD_NON_VALID_BASES_RATIO:
                num_samples += 1
                self.excluded_by_non_valid_bases_ratio += 1
                continue

            sample_gisaid = SampleGisaid(
                run_accession=identifier,
                collection_date=metadata[2],
                host_tax_id=None,
                host=host,
                country_raw=country_raw,
                region=region,
                country=None,
                country_alpha_2=None,
                country_alpha_3=None,
                continent=None,
                continent_alpha_2=None,
                site=None,
                site2=None,
                sequence={CHROMOSOME: compress_sequence(record.seq)},
                sequence_length=sequence_length,
                count_n_bases=count_n_bases,
                count_ambiguous_bases=count_ambiguous_bases,
                covigator_accessor_version=covigator.VERSION,
                status=JobStatus.DOWNLOADED
            )
            try:
                self._parse_country(sample_gisaid)
                self._parse_dates(sample_gisaid)
                samples_gisaid.append(sample_gisaid)
                num_samples += 1
                self.included += 1
                end = time.time()
                total_time += end - start

                if len(samples_gisaid) == BATCH_SIZE:
                    session.add_all(samples_gisaid)
                    session.commit()
                    samples_gisaid = []
            except CovigatorExcludedSampleTooEarlyDateException:
                logger.error("Sample excluded due to too early date")
                self.excluded_by_date += 0

        if len(samples_gisaid) > 0:
            session.add_all(samples_gisaid)
            session.commit()
        if num_samples > 0:
            logger.info("It took {} secs/sample on average".format(float(total_time) / num_samples))
        return num_samples

    def _parse_country(self, gisaid_sample: SampleGisaid):
        parsed_country = self.country_parser.parse_country(gisaid_sample.country_raw)
        gisaid_sample.country = parsed_country.country
        gisaid_sample.country_alpha_2 = parsed_country.country_alpha_2
        gisaid_sample.country_alpha_3 = parsed_country.country_alpha_3
        gisaid_sample.continent_alpha_2 = parsed_country.continent_alpha_2
        gisaid_sample.continent = parsed_country.continent

    def _parse_dates(self, gisaid_sample: SampleGisaid):
        gisaid_sample.collection_date = self._parse_abstract(gisaid_sample.collection_date, date.fromisoformat)
        if gisaid_sample.collection_date is not None and gisaid_sample.collection_date < MINIMUM_DATE:
            raise CovigatorExcludedSampleTooEarlyDateException

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
        logger.info("Excluded by low horizontal coverage = {}".format(self.excluded_by_horizontal_coverage))
        logger.info("Excluded by high non valid bases ratio = {}".format(self.excluded_by_non_valid_bases_ratio))

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
                    "excluded_by_host": self.excluded_by_host,
                    "excluded_by_date": self.excluded_by_date,
                    "excluded_by_horizontal_coverage": self.excluded_by_horizontal_coverage,
                    "excluded_by_non_valid_bases_ratio": self.excluded_by_non_valid_bases_ratio
                }
            }
        ))
        session.commit()
