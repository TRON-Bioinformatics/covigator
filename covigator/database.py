import os
from contextlib import contextmanager

from logzero import logger
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker, Session

from covigator import ENV_COVIGATOR_DB_HOST, ENV_COVIGATOR_DB_NAME, ENV_COVIGATOR_DB_USER, ENV_COVIGATOR_DB_PASSWORD, \
    ENV_COVIGATOR_DB_PORT
from covigator.model import Base


class Database:

    def __init__(self, test=False):
        if test:
            # creates a SQLite in memory database for testing purposes
            db_uri = 'sqlite://'
        else:
            host = os.getenv(ENV_COVIGATOR_DB_HOST, "0.0.0.0")
            database = os.getenv(ENV_COVIGATOR_DB_NAME, "covigator")
            user = os.getenv(ENV_COVIGATOR_DB_USER, "covigator")
            password = os.getenv(ENV_COVIGATOR_DB_PASSWORD, "covigator")
            port = os.getenv(ENV_COVIGATOR_DB_PORT, "5432")
            db_uri = "postgresql+psycopg2://%s:%s@%s:%s/%s" % (user, password, host, port, database)
        self.engine: Engine = create_engine(db_uri)
        self.engine.connect()
        self.Session = sessionmaker(bind=self.engine, autoflush=False)
        self.create_database()

    def create_database(self):
        # this creates all tables in the database (when it exists nothing happens)
        Base.metadata.create_all(self.engine)
        logger.info("Database initialized")

    def get_database_session(self) -> Session:
        return self.Session()


@contextmanager
def session_scope() -> Session:
    """Provide a transactional scope around a series of operations."""
    database = Database()
    session = database.get_database_session()
    try:
        yield session
        session.commit()
    except BaseException as e:
        session.rollback()
        raise e
    finally:
        session.close()
        database.engine.dispose()