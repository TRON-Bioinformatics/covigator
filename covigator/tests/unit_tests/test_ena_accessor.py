from unittest import TestCase
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.database import Database
from covigator.tests import SARS_COV_2_TAXID, HOMO_SAPIENS_TAXID


class FakeEnaAccessor(EnaAccessor):

    def __init__(self, results, database=None):
        # uses an in memory database or the one provided
        super().__init__(tax_id=SARS_COV_2_TAXID, host_tax_id=HOMO_SAPIENS_TAXID,
                         database=database if database else Database(test=True))
        self.results = results

    def _get_ena_runs_page(self, offset):
        return self.results


class EnaAccessorTests(TestCase):

    def test_filtering_by_library_strategies(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "AMPLICON",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "OTHER",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_library_strategy.get("AMPLICON"), 1)
        self.assertEqual(ena_accessor.excluded_samples_by_library_strategy.get("OTHER"), 1)

    def test_filtering_by_instrument_platform(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "OXFORD_NANOPORE",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "SANGER",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_instrument_platform.get("OXFORD_NANOPORE"), 1)
        self.assertEqual(ena_accessor.excluded_samples_by_instrument_platform.get("SANGER"), 1)

    def test_filtering_by_host_taxid(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "1111"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "2222"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_host_tax_id.get("1111"), 1)
        self.assertEqual(ena_accessor.excluded_samples_by_host_tax_id.get("2222"), 1)

    def test_filtering_by_missing_fastqs(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": None,
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_fastq_ftp, 2)

    def test_no_filtering(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)

    def test_filtering_data_already_in_db(self):
        database = Database(test=True)
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ], database=database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)

        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606"},
            {"run_accession": "THIS_IS_A_NEW_ONE",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ], database=database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded_existing, 2)
