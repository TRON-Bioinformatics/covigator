from contextlib import contextmanager
import tenacity
from logzero import logger
from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker, Session
from tenacity import wait_exponential, stop_after_attempt
from covigator.configuration import Configuration
from covigator.database.model import Base, Gene, Conservation
from covigator.exceptions import CovigatorDatabaseConnectionException
from covigator.references.conservation import ConservationLoader
from covigator.references.gene_annotations import GeneAnnotationsLoader


class Database:

    def __init__(self, config: Configuration = None, test=False, verbose=False, initialize=False):
        if test:
            # connects to the test postgres in the CI environment
            #db_uri = "postgresql+psycopg2://%s:%s@%s/%s" % ("test_user", "test_pass", "postgres", "test_db")
            #self.engine: Engine = create_engine(db_uri, echo=verbose)
            initialize = True   # does not make sense not to initialize the DB in the test environment

        db_uri = "postgresql+psycopg2://%s:%s@%s:%s/%s" % (config.db_user, config.db_password, config.db_host,
                                                           config.db_port, config.db_name)
        # these are the default SQLAlchemy values, this values may need to be increased for the dashboard
        self.engine: Engine = create_engine(db_uri, pool_size=config.db_pool_size,
                                            max_overflow=config.db_max_overflow, echo=verbose)
        self.engine.connect()
        self.Session = sessionmaker(bind=self.engine, autoflush=False)
        # avoids initialisation if we know it has already been done
        if initialize:
            self.create_database()

    def create_database(self):
        # this creates all tables in the database (when it exists nothing happens)
        Base.metadata.create_all(self.engine)
        self.initialise_database()
        logger.info("Database initialized")

    def initialise_database(self):
        session = self.get_database_session()
        # loads reference genome if not set
        if session.query(Gene).count() == 0:
            GeneAnnotationsLoader(session).load_data()
        if session.query(Conservation).count() == 0:
            ConservationLoader(session).load_data()

    def get_database_session(self) -> Session:
        return self.Session()


@contextmanager
def session_scope(config: Configuration = None, database: Database = None, initialize=False, test=False) -> Session:
    """Provide a transactional scope around a series of operations.
    Either config or database must be provided!
    """
    if database is None and config is None:
        raise CovigatorDatabaseConnectionException(
            "Must provide either a Configuration object or a database connection to obtain a session")
    if database is None:
        database = Database(config=config, initialize=initialize, test=test)
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


@tenacity.retry(wait=wait_exponential(multiplier=2, min=1, max=10), stop=stop_after_attempt(5))
def get_database(config: Configuration, initialize=False, verbose=False) -> Database:
    try:
        database = Database(config=config, initialize=initialize, verbose=verbose)
        session = database.get_database_session()
        stmt = text("SELECT 1")
        session.execute(stmt)
        logger.info("Database connected!")
    except Exception as e:
        logger.error("Connection to database failed, retrying...")
        raise CovigatorDatabaseConnectionException(e)
    return database
