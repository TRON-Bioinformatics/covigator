from sqlalchemy import and_

from covigator.database.model import SampleEna, JobStatus
from covigator.database import Database
from logzero import logger
from toil.common import Toil
from toil.job import Job
import time


def parent_job(job):
    session = Database().get_database_session()
    count = 0
    try:
        while True:
            ena_run = session.query(SampleEna)\
                .filter(SampleEna.status == JobStatus.PENDING)\
                .order_by(SampleEna.first_created.desc())\
                .first()

            if not ena_run:
                logger.info("No more runs to process after sending {} runs to process".format(count))
                break

            # it has to update the status before doing anything so this processor does not read it again
            ena_run.status = JobStatus.QUEUED
            session.commit()

            # sends the run for processing
            download_job = Job.wrapFn(download, ena_run.run_accession)
            pipeline_job = Job.wrapFn(run_pipeline, ena_run.run_accession)
            delete_job = Job.wrapFn(delete, ena_run.run_accession)
            download_job.addChild(pipeline_job)
            download_job.addFollowOn(delete_job)
            job.addChild(download_job)
            job.addChild(pipeline_job)
    except Exception as e:
        logger.exception(e)
        session.rollback()
    finally:
        session.close()


def download(ena_run_accession: str, memory="2G", cores=1):
    session = Database().get_database_session()
    ena_run = EnaProcessor.find_ena_run_by_accession_and_status(
        ena_run_accession=ena_run_accession, session=session, status=JobStatus.QUEUED)
    if ena_run is not None:
        # TODO: implement download + MD5 check here
        # TODO: we need to configure the location to store the file + store the path in the database
        time.sleep(1)
        logger.info("Downloaded {}: {}".format(ena_run.run_accession, ena_run.fastq_ftp))
        ena_run.status = JobStatus.DOWNLOADED
    else:
        logger.info("Download error {}".format(ena_run_accession))
        ena_run.status = JobStatus.ERROR
    session.commit()


def run_pipeline(ena_run_accession: str, memory="2G", cores=1):
    session = Database().get_database_session()
    ena_run = EnaProcessor.find_ena_run_by_accession_and_status(
        ena_run_accession=ena_run_accession, session=session, status=JobStatus.DOWNLOADED)
    if ena_run is not None:
        # TODO: call the pipeline here
        time.sleep(1)
        logger.info("Processed {}".format(ena_run.run_accession))
        ena_run.status = JobStatus.PROCESSED
    else:
        logger.info("Process error {}".format(ena_run_accession))
        ena_run.status = JobStatus.ERROR
    session.commit()


def delete(ena_run_accession: str, memory="2G", cores=1):
    session = Database().get_database_session()
    ena_run = EnaProcessor.find_ena_run_by_accession_and_status(
        ena_run_accession=ena_run_accession, session=session, status=JobStatus.PROCESSED)
    if ena_run is not None:
        # TODO: delete files from the file system here
        time.sleep(1)
        logger.info("Deleted {}".format(ena_run.run_accession))
        ena_run.status = JobStatus.CLEANED
    else:
        logger.info("Delete error {}".format(ena_run_accession))
        ena_run.status = JobStatus.ERROR
    session.commit()


class EnaProcessor:

    def __init__(self):
        pass

    def process(self, options):
        with Toil(options) as toil:
            toil.start(Job.wrapJobFn(parent_job))

    @staticmethod
    def find_ena_run_by_accession_and_status(ena_run_accession, session, status: JobStatus):
        ena_run = session.query(SampleEna) \
            .filter(and_(SampleEna.run_accession == ena_run_accession, SampleEna.status == status)) \
            .first()
        return ena_run
