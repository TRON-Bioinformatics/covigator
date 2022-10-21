from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from covigator.database.database import Database
from covigator.configuration import Configuration
from covigator.database.model import SampleCovid19Portal



""""
This is a temporary patch that has been amended in the accessor already
"""
config = Configuration()
session = Database(config=config, initialize=True).get_database_session()

sample: SampleCovid19Portal
for sample in session.query(SampleCovid19Portal).all():

    record: SeqRecord
    for record in SeqIO.parse(sample.fasta_path, "fasta"):
        sample.sequence_length = len(str(record.seq))
        sample.count_n_bases = record.seq.count("N")
        sample.count_ambiguous_bases = sum([record.seq.count(b) for b in "RYWSMKHBVD"])
        break

    session.commit()

session.close()
