from datetime import date

from sqlalchemy import and_, desc, asc, Integer
from sqlalchemy.orm import Session

from covigator.database.model import Log, DataSource, CovigatorModule, SampleEna, JobEna, JobStatus


def get_date_of_first_ena_sample(session: Session) -> date:
    """
    Returns the date of the earliest ENA sample loaded in the database
    """
    result = session.query(SampleEna.first_created).join(JobEna).filter(JobEna.status == JobStatus.LOADED) \
        .order_by(asc(SampleEna.first_created)).first()
    return result[0] if result is not None else result


def get_date_of_most_recent_ena_sample(session: Session) -> date:
    """
    Returns the date of the latest ENA sample loaded in the database
    """
    result = session.query(SampleEna.first_created).join(JobEna).filter(JobEna.status == JobStatus.LOADED) \
        .order_by(desc(SampleEna.first_created)).first()
    return result[0] if result is not None else result


def get_date_of_last_check(session: Session, data_source: DataSource) -> date:
    """
    Returns the date of the latest non failed accessor check that also has a subsequent non failed processor run.
    Until the processor has ran the data fetched from the accessor is not available.
    """
    result1 = session.query(Log.start).filter(
        and_(Log.source == data_source, Log.module == CovigatorModule.PROCESSOR, Log.has_error == False)).order_by(
        desc(Log.start)).first()
    most_recent_processor_run = result1[0] if result1 is not None else result1
    result2 = None
    if most_recent_processor_run:
        result2 = session.query(Log.start).filter(
            and_(Log.source == data_source, Log.module == CovigatorModule.ACCESSOR, Log.has_error == False,
                 Log.start < most_recent_processor_run)).order_by(desc(Log.start)).first()
    return result2[0] if result2 is not None else result2


def get_date_of_last_update(session: Session, data_source: DataSource) -> date:
    """
    Returns the date of the latest non failed accessor check **with some new data** that also has a subsequent non
    failed processor run.
    Until the processor has ran the data fetched from the accessor is not available.
    """
    result1 = session.query(Log.start).filter(
        and_(Log.source == data_source, Log.module == CovigatorModule.PROCESSOR, Log.has_error == False)).order_by(
        desc(Log.start)).first()
    most_recent_processor_run = result1[0] if result1 is not None else result1
    result2 = None
    if most_recent_processor_run:
        result2 = session.query(Log.start).filter(
            and_(Log.source == data_source, Log.module == CovigatorModule.ACCESSOR, Log.has_error == False,
                 Log.data["included"].astext.cast(Integer) > 0, Log.start < most_recent_processor_run)) \
            .order_by(desc(Log.start)).first()
    return result2[0] if result2 is not None else result2

