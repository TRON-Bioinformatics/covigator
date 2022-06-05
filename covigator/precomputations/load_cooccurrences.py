from sqlalchemy.orm import Session
from covigator.database.model import DataSource, JobStatus
from covigator.database.queries import Queries
from covigator.pipeline.cooccurrence_matrix import CooccurrenceMatrix
from logzero import logger


class CooccurrenceMatrixLoader:

    def __init__(self, session: Session):
        self.session = session
        self.queries = Queries(session=self.session)
        self.cooccurrence_matrix = CooccurrenceMatrix()

    def load(self, data_source: str, maximum_length: int):

        # deletes the database before loading
        self.session.query(self.queries.get_variant_cooccurrence_klass(data_source)).delete()

        # iterates over every sample in FINISHED status and computes the cooccurrence matrix
        sample_klass = self.queries.get_sample_klass(data_source)
        count_samples = self.queries.count_samples(source=data_source, cache=False)
        computed = 0
        query = self.session.query(sample_klass).filter(sample_klass.status == JobStatus.FINISHED)
        for sample in self.queries.windowed_query(query=query, column=sample_klass.run_accession, windowsize=1000):
            self.cooccurrence_matrix.compute(sample.run_accession, data_source, self.session,
                                             maximum_length=maximum_length)
            computed += 1
            if computed % 1000 == 0:
                logger.info('Processed cooccurrence over {}/{} ({} %) samples'.format(
                    computed, count_samples, round(float(computed) / count_samples * 100, 3)))

        # once finished deletes the unique observations
        variant_cooccurrence_klazz = self.queries.get_variant_cooccurrence_klass(source=data_source)
        self.session.query(self.queries.get_variant_cooccurrence_klass(data_source)) \
            .filter(variant_cooccurrence_klazz.count == 1) \
            .delete()
        self.session.commit()
