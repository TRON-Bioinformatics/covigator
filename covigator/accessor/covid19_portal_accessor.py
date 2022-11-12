import gzip
import os
import pathlib
import shutil
from datetime import datetime, date
from io import StringIO
from json import JSONDecodeError
import random
import time
from typing import Tuple
import pandas as pd

from Bio import SeqIO
from sqlalchemy.exc import IntegrityError

from covigator.configuration import Configuration

from covigator.accessor import MINIMUM_DATE
from requests import Response
from sqlalchemy.orm import Session
import covigator
from covigator.accessor.abstract_accessor import AbstractAccessor, SampleCovid19
from covigator.exceptions import CovigatorExcludedSampleTooEarlyDateException, CovigatorExcludedFailedDownload, \
    CovigatorExcludedTooManyEntries, CovigatorExcludedEmptySequence, \
    CovigatorExcludedBadBases, CovigatorExcludedSampleException, CovigatorExcludedMissingDateException, \
    CovigatorExcludedHorizontalCoverage
from covigator.database.model import DataSource, Log, CovigatorModule, SampleCovid19Portal, JobStatus
from covigator.database.database import Database
from logzero import logger

NUMBER_RETRIES = -1
BATCH_SIZE = 1000
THRESHOLD_NON_VALID_BASES_RATIO = 0.2
THRESHOLD_GENOME_COVERAGE = 0.2
GENOME_LENGTH = 29903


class Covid19PortalAccessor(AbstractAccessor):

    API_URL_BASE = "https://www.covid19dataportal.org/api/backend/viral-sequences/sequences"
    FASTA_URL_BASE = "https://www.ebi.ac.uk/ena/browser/api/fasta"
    # while ENA API is in principle organism agnostic, this is SARS-CoV-2 specific, thus taxon is hard-coded
    HOST = "Homo sapiens"
    HOST2 = "Human"
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
        self.excluded_missing_date = 0
        self.excluded_failed_download = 0
        self.excluded_empty_sequence = 0
        self.excluded_too_many_entries = 0
        self.excluded_horizontal_coverage = 0
        self.excluded_bad_bases = 0

    def access(self):
        logger.info("Starting Covid19 portal accessor")
        session = self.database.get_database_session()

        # NOTE: holding in memory the whole list of existing ids is much faster than querying every time
        # using a set is way faster than a list
        existing_sample_ids = set([value for value, in session.query(SampleCovid19Portal.run_accession).all()])
        try:
            logger.info("Reading...")
            count = 0

            # queries data by month since December 2019, the API does not support paginating over more than 1M entries
            # there is no month with more than 1M entries... at least for now...
            months = pd.date_range('2019-12-01', date.today().strftime("%Y-%m-%d"), freq='MS').strftime("%Y%m").tolist()
            for m in months:
                page = 0
                logger.info("Querying for month {}".format(m))
                while True:
                    page += 1   # page starts at 1
                    status_code, list_runs = self._get_page(page=page, size=BATCH_SIZE, month=m)
                    entries = list_runs.get('entries')
                    if entries is None or entries == []:
                        logger.info("Last page reached!")
                        logger.info("Status code: {}".format(status_code))
                        logger.info("Content: {}".format(str(list_runs)))
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

    def _get_page(self, page, size, month) -> Tuple[int, dict]:
        # the API is sometimes unstable returning bad json, we retry
        success = False
        json = None
        status_code = None
        retries_count = 0
        backoff_iteration = 1
        truncate_iteration = 8
        while not success:
            response: Response = self.get_with_retries(
                "{url_base}?page={page}&size={size}&query=collection_date:({month})&fields=lineage,coverage,collection_date,country,host,TAXON,"
                "creation_date,last_modification_date,center_name,isolate,molecule_type".format(
                    url_base=self.API_URL_BASE, page=page, size=size, month=month))
            try:
                json = response.json()
                status_code = response.status_code
                success = True
            except JSONDecodeError:
                # retry
                retries_count += 1
                # waits for an increasing random time
                random_sleep = random.randrange(0, (2 ** backoff_iteration) - 1)
                logger.info("Retrying connection after %s seconds" % str(random_sleep))
                time.sleep(random_sleep)
                # when it reaches the maximum value that it may wait it stops increasing time
                if backoff_iteration < truncate_iteration:
                    backoff_iteration += 1
        return status_code, json

    def _process_runs(self, list_samples, existing_sample_ids, session: Session):

        included_samples = []
        for sample_dict in list_samples.get('entries'):
            if isinstance(sample_dict, dict):
                run_accession = self._get_run_accession(sample_dict)
                if run_accession in existing_sample_ids:
                    self.excluded_existing += 1
                    continue    # skips runs already registered in the database
                if not self._complies_with_inclusion_criteria(sample_dict):
                    continue    # skips runs not complying with inclusion criteria
                # NOTE: this parse operation is costly
                try:
                    # parses sample into DB model
                    sample = self._parse_covid19_portal_sample(sample_dict)
                    # downloads FASTA file
                    sample = self._download_fasta(sample=sample)
                    self.included += 1
                    included_samples.append(sample)
                except CovigatorExcludedSampleTooEarlyDateException:
                    self.excluded_by_date += 1
                    self.excluded += 1
                except CovigatorExcludedMissingDateException:
                    self.excluded_missing_date += 1
                    self.excluded += 1
                except CovigatorExcludedFailedDownload:
                    self.excluded_failed_download += 1
                    self.excluded += 1
                except CovigatorExcludedEmptySequence:
                    self.excluded_empty_sequence += 1
                    self.excluded += 1
                except CovigatorExcludedTooManyEntries:
                    self.excluded_too_many_entries += 1
                    self.excluded += 1
                except CovigatorExcludedHorizontalCoverage:
                    self.excluded_horizontal_coverage += 1
                    self.excluded += 1
                except CovigatorExcludedBadBases:
                    self.excluded_bad_bases += 1
                    self.excluded += 1
                except CovigatorExcludedSampleException:
                    self.excluded += 1
            else:
                logger.error("Sample without the expected format")

        if len(included_samples) > 0:
            try:
                session.add_all(included_samples)
                session.commit()
            except IntegrityError:
                # NOTE: we observed that the API may return repeated samples, thus when this error is raised,
                # we insert one by one and ignore the failing sample
                # this is slower but it happens rarely
                session.rollback()
                for s in included_samples:
                    try:
                        session.add(s)
                        session.commit()
                    except IntegrityError:
                        pass
        self._log_results()

    def _parse_covid19_portal_sample(self, sample: dict) -> SampleCovid19Portal:
        run_accession = self._get_run_accession(sample)
        if run_accession is None:
            raise CovigatorExcludedSampleException("Missing sample id")
        sample = SampleCovid19Portal(
            run_accession=run_accession,  # this used to be "acc"
            first_created=next(iter(sample.get('fields').get('creation_date')), None),
            collection_date=next(iter(sample.get('fields').get('collection_date')), None),
            last_modification_date=next(iter(sample.get('fields').get('last_modification_date')), None),
            center_name=next(iter(sample.get('fields').get('center_name')), None),
            isolate=next(iter(sample.get('fields').get('isolate')), None),
            molecule_type=next(iter(sample.get('fields').get('molecule_type')), None),
            country=next(iter(sample.get('fields').get('country')), None),
            pangolin_lineage=next(iter(sample.get('fields').get('lineage')), None),
            # build FASTA URL here
            fasta_url="{base}/{acc}".format(base=self.FASTA_URL_BASE, acc=run_accession)
        )
        self._parse_country(sample)
        self._parse_dates(sample)
        sample.covigator_accessor_version = covigator.VERSION
        return sample

    def _get_run_accession(self, sample):
        run_accession = sample.get('id')
        if run_accession is None:
            run_accession = sample.get('acc')
        return run_accession

    def _complies_with_inclusion_criteria(self, sample: dict):
        # NOTE: this uses the original dictionary instead of the parsed SampleEna class for performance reasons
        included = True

        host = next(iter(sample.get('fields').get('host')), None)
        if host is None or not self._match_host(host):
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

    def _match_host(self, host: str):
        return host.startswith(self.HOST) or host.startswith(self.HOST2)

    def _log_results(self):
        logger.info("Included new runs = {}".format(self.included))
        logger.info("Excluded already existing samples = {}".format(self.excluded_existing))
        logger.info("Total excluded runs by selection criteria = {}".format(self.excluded))
        logger.info("Excluded by host = {}".format(self.excluded_samples_by_host_tax_id))
        logger.info("Excluded by host = {}".format(self.excluded_samples_by_tax_id))
        logger.info("Excluded failed download = {}".format(self.excluded_failed_download))
        logger.info("Excluded empty sequence = {}".format(self.excluded_empty_sequence))
        logger.info("Excluded horizontal coverage = {}".format(self.excluded_horizontal_coverage))
        logger.info("Excluded too many entries = {}".format(self.excluded_too_many_entries))
        logger.info("Excluded bad bases = {}".format(self.excluded_bad_bases))
        logger.info("Excluded missing date = {}".format(self.excluded_missing_date))

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
                    "excluded_early_date": self.excluded_by_date,
                    "excluded_missing_date": self.excluded_missing_date,
                    "excluded_failed_download": self.excluded_failed_download,
                    "excluded_bad_bases": self.excluded_bad_bases,
                    "excluded_horizontal_coverage": self.excluded_horizontal_coverage,
                    "excluded_too_many_entries": self.excluded_too_many_entries,
                    "excluded_empty_sequence": self.excluded_empty_sequence,
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

        # avoids downloading the same files over and over
        if not os.path.exists(local_full_path):
            pathlib.Path(local_folder).mkdir(parents=True, exist_ok=True)
            try:
                fasta_str = self.get_with_retries(sample.fasta_url).content.decode('utf-8')
                fasta_io = StringIO(fasta_str)
                records = list(SeqIO.parse(fasta_io, "fasta"))
            except Exception as e:
                raise CovigatorExcludedFailedDownload(e)
        else:
            records = list(SeqIO.parse(open(local_full_path, "r"), "fasta"))

        # checks the validity of the FASTA sequence
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

        # compress and writes the file after all checks
        if not os.path.exists(local_full_path):
            with open(local_full_path, "w") as f:
                SeqIO.write(sequences=records, handle=f, format="fasta")

        # stores the reference to the file in the DB and metadata about the sequence
        sample.fasta_path = local_full_path
        sample.sequence_length = sequence_length
        sample.count_n_bases = count_n_bases
        sample.count_ambiguous_bases = count_ambiguous_bases
        sample.status = JobStatus.DOWNLOADED

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
        except Exception:
            sample.collection_date = None
        try:
            sample.first_created = datetime.strptime(sample.first_created, format).date()
        except Exception:
            sample.first_created = None

        if sample.collection_date is None:
            raise CovigatorExcludedMissingDateException()
        if sample.collection_date is not None and sample.collection_date < MINIMUM_DATE:
            raise CovigatorExcludedSampleTooEarlyDateException()


if __name__ == '__main__':
    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_accesor)
    Covid19PortalAccessor(database=Database(config=config, initialize=True), storage_folder=config.storage_folder).access()
