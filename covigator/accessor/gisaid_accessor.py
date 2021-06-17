import base64
from datetime import date, datetime
import pycountry
import pycountry_convert
from sqlalchemy.orm import Session
from covigator.database.model import SampleGisaid, JobGisaid, Sample, DataSource, Log, CovigatorModule
from covigator.database.database import Database
from logzero import logger
from Bio import SeqIO
import time
import lz4.frame

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
        self.excluded_existing = 0
        self.included = 0
        self.cache_countries = {}

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
        jobs_and_samples = []
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

            compressed_sequence = self.compress_sequence(record.seq)
            sample_gisaid = SampleGisaid(
                run_accession=identifier,
                date=metadata[2],
                host_tax_id=None,
                host=host,
                country_raw=country_raw,
                region=region,
                country=None,
                country_alpha_2=None,
                country_alpha_3=None,
                continent=location_list[0],
                continent_alpha_2=None,
                site=None,
                site2=None,
                sequence={"MN908947.3": compressed_sequence}
            )

            self._parse_country(sample_gisaid)
            self._parse_dates(sample_gisaid)
            sample = self._build_sample(sample_gisaid)
            job = JobGisaid(run_accession=sample_gisaid.run_accession)
            samples_gisaid.append(sample_gisaid)
            jobs_and_samples.append(job)
            jobs_and_samples.append(sample)
            num_samples += 1
            end = time.time()
            total_time += end - start

            if len(samples_gisaid) == BATCH_SIZE:
                session.add_all(samples_gisaid)
                session.commit()
                session.add_all(jobs_and_samples)
                session.commit()
                samples_gisaid = []
                jobs_and_samples = []

        if len(samples_gisaid) > 0:
            session.add_all(samples_gisaid)
            session.commit()
            session.add_all(jobs_and_samples)
            session.commit()
        logger.info("It took {} secs/sample on average".format(float(total_time) / num_samples))
        return num_samples

    @staticmethod
    def compress_sequence(sequence):
        return base64.b64encode(lz4.frame.compress(sequence.encode('utf-8'))).decode()

    def _build_sample(self, sample_gisaid):
        return Sample(
            id=sample_gisaid.run_accession,
            source=DataSource.GISAID,
            gisaid_id=sample_gisaid.run_accession
        )

    def _parse_country(self, gisaid_sample: SampleGisaid):
        try:
            if gisaid_sample.country_raw is not None and gisaid_sample.country_raw.strip() == "":
                gisaid_sample.country_raw = None
                raise LookupError   # if no input avoids trying to parse it
            match = self.cache_countries.get(gisaid_sample.country_raw)
            if match is None:
                match = pycountry.countries.search_fuzzy(gisaid_sample.country_raw)[0]
                self.cache_countries[gisaid_sample.country_raw] = match
            gisaid_sample.country = match.name
            gisaid_sample.country_alpha_2 = match.alpha_2
            gisaid_sample.country_alpha_3 = match.alpha_3
            gisaid_sample.continent_alpha_2 = pycountry_convert.country_alpha2_to_continent_code(gisaid_sample.country_alpha_2)
            gisaid_sample.continent = pycountry_convert.convert_continent_code_to_continent_name(gisaid_sample.continent_alpha_2)
        except (LookupError, AttributeError):
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
