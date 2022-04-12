import traceback
from contextlib import suppress
from datetime import datetime
from logzero import logger
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.database.model import Log, DataSource, CovigatorModule, JobStatus
from covigator.database.queries import Queries
from covigator.misc import backoff_retrier
from covigator.pipeline.downloader import Downloader


class EnaDownloader:

    def __init__(self, database: Database, config: Configuration, wait_time=60):
        self.config = config
        self.start_time = datetime.now()
        self.has_error = False
        self.error_message = None
        self.database = database
        assert self.database is not None, "Empty database"
        self.session = self.database.get_database_session()
        self.queries = Queries(self.session)
        self.wait_time = wait_time
        self.downloader = Downloader(config=config)

    def process(self):
        logger.info("Starting ENA downloader")
        count = 0
        try:
            while True:
                # queries 100 jobs every time to make sending to queue faster
                samples = self.queries.find_first_jobs_to_download(DataSource.ENA, n=1000)
                if samples is None or len(samples) == 0:
                    logger.info("No more ENA samples to download after downloading {} samples".format(count))
                    break
                for sample in samples:
                    # downloads
                    try:
                        download_with_retries = backoff_retrier.wrapper(self.downloader.download, self.config.retries_download)
                        fastq_path = download_with_retries(sample_ena=sample)
                        sample.fastq_path = fastq_path
                        sample.downloaded_at = datetime.now()
                        sample.sample_folder = sample.get_sample_folder(base_folder=self.config.storage_folder)
                        sample.status = JobStatus.DOWNLOADED
                        self.session.commit()
                    except Exception as e:
                        sample.status = JobStatus.FAILED_DOWNLOAD
                        sample.error_message = self._get_traceback_from_exception(e)
                        self.session.commit()

                    count += 1
                    if count % 1000 == 0:
                        logger.info("Downloaded {} ENA samples...".format(count))

            logger.info("ENA downloader finished!")

        except Exception as e:
            logger.exception(e)
            self.session.rollback()
            self.has_error = True
        finally:
            logger.info("Logging execution stats...")
            self._write_execution_log(count, data_source=DataSource.ENA)
            with suppress(Exception):
                self.session.close()
            logger.info("Database sessions closed")

    def _write_execution_log(self, count, data_source: DataSource):
        end_time = datetime.now()
        self.session.add(Log(
            start=self.start_time,
            end=end_time,
            source=data_source,
            module=CovigatorModule.DOWNLOADER,
            has_error=self.has_error,
            processed=count,
            data={
                "processed": count
            }
        ))
        self.session.commit()

    @staticmethod
    def _get_traceback_from_exception(e):
        return "".join(traceback.format_exception(
            etype=type(e), value=e, tb=e.__traceback__))

