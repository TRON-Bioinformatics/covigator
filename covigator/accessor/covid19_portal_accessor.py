import gzip
import os
import pathlib
import shutil
from datetime import datetime
from json import JSONDecodeError
from urllib import request
from Bio import SeqIO
from covigator.accessor import MINIMUM_DATE

from covigator.configuration import Configuration
from requests import Response
from sqlalchemy.orm import Session
import covigator
from covigator.accessor.abstract_accessor import AbstractAccessor, SampleCovid19
from covigator.exceptions import CovigatorExcludedSampleTooEarlyDateException, CovigatorExcludedFailedDownload, \
    CovigatorExcludedTooManyEntries, CovigatorExcludedEmptySequence, CovigatorExcludedHorizontalCoverage, \
    CovigatorExcludedBadBases, CovigatorExcludedSampleException, CovigatorExcludedMissingDateException
from covigator.database.model import DataSource, Log, CovigatorModule, SampleCovid19Portal
from covigator.database.database import Database
from logzero import logger

from covigator.misc import backoff_retrier

NUMBER_RETRIES = 5
BATCH_SIZE = 1000
THRESHOLD_NON_VALID_BASES_RATIO = 0.2
THRESHOLD_GENOME_COVERAGE = 0.2
GENOME_LENGTH = 29903


class Covid19PortalAccessor(AbstractAccessor):

    API_URL_BASE = "https://www.covid19dataportal.org/api/backend/viral-sequences/sequences"
    FASTA_URL_BASE = "https://www.ebi.ac.uk/ena/browser/api/fasta"
    # while ENA API is in principle organism agnostic, this is SARS-CoV-2 specific, thus taxon is hard-coded
    HOST = "Homo sapiens"
    TAX_ID = "2697049"

    def __init__(self, database: Database, storage_folder):
        super().__init__()
        logger.info("Initialising Covid19 portal accessor")
        self.start_time = datetime.now()
        self.has_error = False
        self.error_message = None
        self.storage_folder = storage_folder

        self.database = database
        assert self.database is not None, "Empty database"

        self.excluded_samples_by_host_tax_id = {}
        self.excluded_samples_by_tax_id = {}
        self.excluded_existing = 0
        self.included = 0
        self.excluded = 0
        self.excluded_by_date = 0
        self.excluded_failed_download = 0

        self.download_with_retries = backoff_retrier.wrapper(self._download_fasta, NUMBER_RETRIES)

    def access(self):
        logger.info("Starting Covid19 portal accessor")
        session = self.database.get_database_session()
        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # it assumes there will be no repetitions
        existing_sample_ids = [value for value, in session.query(SampleCovid19Portal.run_accession).all()]
        try:
            logger.info("Reading...")
            page = 1    # it has to start with 1
            list_runs = self._get_page(page=page, size=BATCH_SIZE)
            num_entries = len(list_runs.get('entries'))
            count = num_entries
            # gets total expected number of results from first page
            logger.info("Processing {} Covid19 Portal samples...".format(num_entries))
            self._process_runs(list_runs, existing_sample_ids, session)

            while True:
                page += 1
                list_runs = self._get_page(page=page, size=BATCH_SIZE)
                entries = list_runs.get('entries')
                if entries is None or entries == []:
                    break
                num_entries = len(entries)
                count += num_entries
                self._process_runs(list_runs, existing_sample_ids, session)
                logger.info("Processed {} of Covid19 Portal samples...".format(count))

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
            logger.info("Finished Covid19 Portal accessor")

    def _get_page(self, page, size) -> dict:
        # as communicated by ENA support we use limit=0 and offset=0 to get all records in one query
        response: Response = self.get_with_retries(
            "{url_base}?page={page}&size={size}".format(url_base=self.API_URL_BASE, page=page, size=size))
        try:
            json = response.json()
        except JSONDecodeError as e:
            logger.exception(e)
            raise e
        return json

    def _process_runs(self, list_samples, existing_sample_ids, session: Session):

        included_samples = []
        for sample_dict in list_samples.get('entries'):
            if isinstance(sample_dict, dict):
                if sample_dict.get("acc") in existing_sample_ids:
                    self.excluded_existing += 1
                    continue    # skips runs already registered in the database
                if not self._complies_with_inclusion_criteria(sample_dict):
                    continue    # skips runs not complying with inclusion criteria
                # NOTE: this parse operation is costly
                try:
                    # parses sample into DB model
                    sample = self._parse_covid19_portal_sample(sample_dict)
                    # downloads FASTA file
                    sample = self.download_with_retries(sample=sample)
                    self.included += 1
                    included_samples.append(sample)
                except CovigatorExcludedSampleTooEarlyDateException:
                    self.excluded_by_date += 1
                    self.excluded += 1
                except CovigatorExcludedFailedDownload:
                    self.excluded_failed_download += 1
                    self.excluded += 1
                except CovigatorExcludedSampleException:
                    self.excluded += 1
            else:
                logger.error("Sample without the expected format")
                logger.error("Sample without the expected format")

        if len(included_samples) > 0:
            session.add_all(included_samples)
            session.commit()

    def _parse_covid19_portal_sample(self, sample: dict) -> SampleCovid19Portal:
        sample = SampleCovid19Portal(
            run_accession=sample.get('acc'),
            first_created=next(iter(sample.get('fields').get('creation_date')), None),
            collection_date=next(iter(sample.get('fields').get('collection_date')), None),
            last_modification_date=next(iter(sample.get('fields').get('last_modification_date')), None),
            center_name=next(iter(sample.get('fields').get('center_name')), None),
            isolate=next(iter(sample.get('fields').get('isolate')), None),
            molecule_type=next(iter(sample.get('fields').get('molecule_type')), None),
            country=next(iter(sample.get('fields').get('country')), None),
            # build FASTA URL here
            fasta_url="{base}/{acc}".format(base=self.FASTA_URL_BASE, acc=sample.get('acc'))
        )
        self._parse_country(sample)
        self._parse_dates(sample)
        sample.covigator_accessor_version = covigator.VERSION
        return sample

    def _complies_with_inclusion_criteria(self, sample: dict):
        # NOTE: this uses the original dictionary instead of the parsed SampleEna class for performance reasons
        included = True

        host = next(iter(sample.get('fields').get('host')), None)
        if host is None or host.strip() == "" or host != self.HOST:
            included = False    # skips runs where the host is empty or does not match
            self.excluded_samples_by_host_tax_id[str(host)] = \
                self.excluded_samples_by_host_tax_id.get(str(host), 0) + 1

        taxon = next(iter(sample.get('fields').get('TAXON')), None)
        if taxon is None or taxon.strip() == "" or taxon != self.TAX_ID:
            included = False  # skips runs where the host is empty or does not match
            self.excluded_samples_by_tax_id[str(taxon)] = \
                self.excluded_samples_by_tax_id.get(str(taxon), 0) + 1

        if not included:
            self.excluded += 1
        return included

    def _log_results(self):
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded already existing samples = {}".format(self.excluded_existing))
        logger.info("Total excluded runs by selection criteria = {}".format(self.excluded))
        logger.info("Excluded by host = {}".format(self.excluded_samples_by_host_tax_id))
        logger.info("Excluded by host = {}".format(self.excluded_samples_by_tax_id))
        logger.info("Excluded failed download = {}".format(self.excluded_failed_download))

    def _write_execution_log(self, session: Session):
        end_time = datetime.now()
        session.add(Log(
            start=self.start_time,
            end=end_time,
            source=DataSource.COVID19_PORTAL,
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
                    "excluded_failed_download": self.excluded_failed_download,
                    "host": self.excluded_samples_by_host_tax_id,
                    "taxon": self.excluded_samples_by_tax_id
                }
            }
        ))
        session.commit()

    def _download_fasta(self, sample: SampleCovid19Portal) -> SampleCovid19Portal:

        local_filename = sample.fasta_url.split('/')[-1] + ".fasta"  # URL comes without extension
        local_folder = sample.get_sample_folder(self.storage_folder)
        local_full_path = os.path.join(local_folder, local_filename)
        compressed_local_full_path = local_full_path + ".gz"

        # avoids downloading the same files over and over
        if not os.path.exists(local_full_path):
            pathlib.Path(local_folder).mkdir(parents=True, exist_ok=True)
            try:
                request.urlretrieve(sample.fasta_url, local_full_path)
                self._compress_file(local_full_path, compressed_local_full_path)
            except Exception as e:
                raise CovigatorExcludedFailedDownload(e)

        sample.fasta_path = compressed_local_full_path

        # checks the validity of the FASTA sequence
        records = list(SeqIO.parse(gzip.open(compressed_local_full_path, "rt"), "fasta"))
        if len(records) == 0:
            raise CovigatorExcludedEmptySequence()
        if len(records) > 1:
            raise CovigatorExcludedTooManyEntries()
        record = records[0]
        sequence_length = len(record.seq)
        if float(sequence_length) / GENOME_LENGTH < THRESHOLD_GENOME_COVERAGE:
            raise CovigatorExcludedHorizontalCoverage()
        count_n_bases = record.seq.count("N")
        count_ambiguous_bases = sum([record.seq.count(b) for b in "RYWSMKHBVD"])
        if float(count_n_bases + count_ambiguous_bases) / sequence_length > THRESHOLD_NON_VALID_BASES_RATIO:
            raise CovigatorExcludedBadBases()

        return sample

    def _compress_file(self, uncompressed_file, compressed_file):
        with open(uncompressed_file, 'rb') as f_in:
            with gzip.open(compressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(uncompressed_file)

    def _parse_dates(self, sample: SampleCovid19):
        format = "%Y%m%d"
        try:
            sample.collection_date = datetime.strptime(sample.collection_date, format).date()
        except ValueError:
            sample.collection_date = None
        try:
            sample.first_created = datetime.strptime(sample.first_created, format).date()
        except ValueError:
            sample.first_created = None

        if sample.collection_date is None:
            raise CovigatorExcludedMissingDateException()
        if sample.collection_date is not None and sample.collection_date < MINIMUM_DATE:
            raise CovigatorExcludedSampleTooEarlyDateException()


if __name__ == '__main__':
    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_accesor)
    Covid19PortalAccessor(database=Database(config=config, initialize=True), storage_folder=config.storage_folder).access()