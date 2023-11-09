from datetime import date
from covigator.database.model import SampleEna, Log, DataSource, CovigatorModule
from covigator.tests.unit_tests.abstract_test import AbstractTest
from covigator.tests.unit_tests.faked_objects import FakeEnaAccessor


class EnaAccessorTests(AbstractTest):

    def test_filtering_by_library_strategies(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "AMPLICON",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "OTHER",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
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
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "SANGER",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
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
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "1111"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "2222"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_host_tax_id.get("1111"), 1)
        self.assertEqual(ena_accessor.excluded_samples_by_host_tax_id.get("2222"), 1)

    def test_filtering_by_host_taxid_disabled(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080484_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "1111"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "2222"},
            {"run_accession": "ERR4080486",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "OTHER",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/004/ERR4080484/ERR4080486_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "OXFORD_NANOPORE",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"}
        ], host_tax_id=None)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_instrument_platform.get("OXFORD_NANOPORE"), 1)
        self.assertEqual(ena_accessor.excluded_samples_by_library_strategy.get("OTHER"), 1)

    def test_filtering_by_minimum_date(self):
        ena_accessor_with_min_date = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
            "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
            "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "first_created": "2020-01-01",
             "collection_date": "2019-12-31",
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "first_created": "2020-01-01 14:50",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "first_created": "2018-12-31",
             "collection_date": "2018-12-30"}
        ], disable_minimum_date=False)
        ena_accessor_with_min_date.access()
        self.assertEqual(ena_accessor_with_min_date.included, 2)
        self.assertEqual(ena_accessor_with_min_date.excluded_by_date, 1)

        ena_accessor_without_min_date = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080486",
            "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
            "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080486/ERR4080486_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "first_created": "2020-01-01",
             "collection_date": "2019-12-31",
             },
            {"run_accession": "ERR4080487",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080487/ERR4080487_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "first_created": "2020-01-01 14:50",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080488",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080489/ERR4080489_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "first_created": "2018-12-31",
             "collection_date": "2018-12-30"}
        ], disable_minimum_date=True)
        ena_accessor_without_min_date.access()
        self.assertEqual(ena_accessor_without_min_date.included, 3)
        self.assertEqual(ena_accessor_without_min_date.excluded_by_date, 0)

    def test_filtering_by_missing_fastqs(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": None,
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded, 2)
        self.assertEqual(ena_accessor.excluded_samples_by_fastq_ftp, 2)

    def test_filtering_empty_collection_date(self):
        # Test that samples with empty collection date are filtered out
        ena_accessor_empty_collection = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
            "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
            "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "first_created": "2020-01-01",
             "collection_date": "2019-12-31",
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "first_created": "2020-01-01 14:50",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"
            }
        ], disable_collection_date=False)
        ena_accessor_empty_collection.access()
        self.assertEqual(ena_accessor_empty_collection.included, 2)
        self.assertEqual(ena_accessor_empty_collection.excluded, 1)
        self.assertEqual(ena_accessor_empty_collection.excluded_samples_by_empty_collection_date, 1)

        # Test that collection date filter can be disabled
        ena_accessor_empty_collection_disabled = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080486",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080486_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "first_created": "2020-01-01",
             "collection_date": "2019-12-31 12:12:12",
             },
            {"run_accession": "ERR4080487",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080487_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "first_created": "2020-01-01 14:50",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080488",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080488_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "first_created": "",
             "collection_date": ""}
        ], disable_collection_date=True)
        ena_accessor_empty_collection_disabled.access()
        self.assertEqual(ena_accessor_empty_collection_disabled.included, 3)
        self.assertEqual(ena_accessor_empty_collection_disabled.exlcuded, 0)
        self.assertEqual(ena_accessor_empty_collection_disabled.excluded_samples_by_empty_collection_date, 0)

    def test_no_filtering(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)

    def test_filtering_data_already_in_db(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "collection_date": "2019-12-31 12:12:12",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606"}
        ], database=self.database)
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
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "THIS_IS_A_NEW_ONE",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"}
        ], database=self.database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 1)
        self.assertEqual(ena_accessor.excluded_existing, 2)

    def test_country_parsing(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12",
             "country": "england"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12",
             "country": "GermaN"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12",
             "country": "Morocco:Meknez"},
            {"run_accession": "ERR4080486",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "collection_date": "2019-12-31 12:12:12",
             "host_tax_id": "9606"},
            {"run_accession": "ERR4080487",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12",
             "country": ""},
            {"run_accession": "ERR4080488",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12",
             "country": "Jupiter"}
        ], database=self.database, disable_collection_date=True)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 6)
        self.assertEqual(ena_accessor.excluded, 0)
        session = self.database.get_database_session()
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080483").first()
        self.assertEqual(run.country_raw, "england")
        self.assertEqual(run.country, "United Kingdom")
        self.assertEqual(run.country_alpha_2, "GB")
        self.assertEqual(run.country_alpha_3, "GBR")
        self.assertEqual(run.continent_alpha_2, "EU")
        self.assertEqual(run.continent, "Europe")
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080484").first()
        self.assertEqual(run.country_raw, "GermaN")
        self.assertEqual(run.country, "Germany")
        self.assertEqual(run.country_alpha_2, "DE")
        self.assertEqual(run.country_alpha_3, "DEU")
        self.assertEqual(run.continent_alpha_2, "EU")
        self.assertEqual(run.continent, "Europe")
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080485").first()
        self.assertEqual(run.country_raw, "Morocco:Meknez")
        self.assertEqual(run.country, "Morocco")
        self.assertEqual(run.country_alpha_2, "MA")
        self.assertEqual(run.country_alpha_3, "MAR")
        self.assertEqual(run.continent_alpha_2, "AF")
        self.assertEqual(run.continent, "Africa")
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080486").first()
        self.assertEqual(run.country_raw, None)
        self.assertEqual(run.country, "Not available")
        self.assertEqual(run.country_alpha_2, "None")
        self.assertEqual(run.country_alpha_3, "None")
        self.assertEqual(run.continent_alpha_2, "None")
        self.assertEqual(run.continent, "None")
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080487").first()
        self.assertEqual(run.country_raw, "")
        self.assertEqual(run.country, "Not available")
        self.assertEqual(run.country_alpha_2, "None")
        self.assertEqual(run.country_alpha_3, "None")
        self.assertEqual(run.continent_alpha_2, "None")
        self.assertEqual(run.continent, "None")
        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080488").first()
        self.assertEqual(run.country_raw, "Jupiter")
        self.assertEqual(run.country, "Not available")
        self.assertEqual(run.country_alpha_2, "None")
        self.assertEqual(run.country_alpha_3, "None")
        self.assertEqual(run.continent_alpha_2, "None")
        self.assertEqual(run.continent, "None")

    def test_dates_parsing(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "first_created": "2020-01-01",
             "collection_date": "2019-12-31",
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "first_created": "2020-01-01 14:50",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "first_created": "blah",
             "collection_date": "blah"}
        ], database=self.database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)
        session = self.database.get_database_session()

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080483").first()
        self.assertEqual(run.collection_date, date.fromisoformat("2019-12-31"))
        self.assertEqual(run.first_created, date.fromisoformat("2020-01-01"))

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080484").first()
        self.assertIsNone(run.collection_date)
        self.assertIsNone(run.first_created)

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080485").first()
        self.assertIsNone(run.collection_date)
        self.assertIsNone(run.first_created)

    def test_numeric_values(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "lat": "",
             "lon": "",
             "read_count": "",
             "base_count": "",
             "nominal_length": "",
             "collection_date": "2019-12-31 12:12:12"
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "lat": "hey",
             "lon": "hey",
             "read_count": "hey",
             "base_count": "hey",
             "nominal_length": "hey",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "lat": "1.1",
             "lon": "1.1",
             "read_count": "1",
             "base_count": "1",
             "nominal_length": "1",
             "collection_date": "2019-12-31 12:12:12"
             }
        ], database=self.database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)
        session = self.database.get_database_session()

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080483").first()
        self.assertIsNone(run.lat)
        self.assertIsNone(run.lon)
        self.assertIsNone(run.nominal_length)
        self.assertIsNone(run.read_count)
        self.assertIsNone(run.base_count)

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080484").first()
        self.assertIsNone(run.lat)
        self.assertIsNone(run.lon)
        self.assertIsNone(run.nominal_length)
        self.assertIsNone(run.read_count)
        self.assertIsNone(run.base_count)

        run = session.query(SampleEna).filter(SampleEna.run_accession == "ERR4080485").first()
        self.assertEqual(run.lat, 1.1)
        self.assertEqual(run.lon, 1.1)
        self.assertEqual(run.nominal_length, 1)
        self.assertEqual(run.read_count, 1)
        self.assertEqual(run.base_count, 1)

    def test_sample_and_job_loading(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
            },
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
             }
        ], database=self.database)
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 3)
        self.assertEqual(ena_accessor.excluded, 0)
        session = self.database.get_database_session()

        self._assert_entities_from_accessor(session, "ERR4080483")
        self._assert_entities_from_accessor(session, "ERR4080484")
        self._assert_entities_from_accessor(session, "ERR4080485")

    def _assert_entities_from_accessor(self, session, identifier):
        self.assertEqual(session.query(SampleEna).filter(SampleEna.run_accession == identifier).count(), 1)

    def test_writing_logs(self):
        ena_accessor = FakeEnaAccessor(results=[
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
             },
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/003/ERR4080483/ERR4080483_1.fastq.gz",
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
            },
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"
             }
        ], database=self.database)
        ena_accessor.access()

        session = self.database.get_database_session()
        self.assertEqual(session.query(Log).count(), 1)
        log = session.query(Log).first()
        self.assertIsNotNone(log.start)
        self.assertIsNotNone(log.end)
        self.assertEqual(log.source, DataSource.ENA)
        self.assertEqual(log.module, CovigatorModule.ACCESSOR)
        self.assertEqual(log.processed, 3)
        data = log.data
        self.assertEqual(data.get("included"), 3)
        self.assertEqual(data.get("excluded").get("existing"), 0)
        self.assertEqual(data.get("excluded").get("excluded_by_criteria"), 0)
        self.assertEqual(data.get("excluded").get("missing_fastq"), 0)
        self.assertIsInstance(data.get("excluded").get("platform"), dict)
        self.assertIsInstance(data.get("excluded").get("library_strategy"), dict)
        self.assertIsInstance(data.get("excluded").get("host"), dict)

    def test_excluding_rnaseq_samples(self):
        ena_accessor = FakeEnaAccessor([
            {"run_accession": "ERR4080483",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "ILLUMINA",
             "library_strategy": "RNA-Seq",
             "fastq_ftp": "",
             "fastq_md5": "a91a9dfa2f7008e13a7ce9767aa9aaf3",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080484",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "RNA-Seq",
             "library_strategy": "WGS",
             "fastq_ftp": None,
             "fastq_md5": "c57fef34933cbbec2e9e08867f3c664c",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"},
            {"run_accession": "ERR4080485",
             "scientific_name": "Severe acute respiratory syndrome coronavirus 2",
             "instrument_platform": "RNA-Seq",
             "library_strategy": "WGS",
             "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR408/005/ERR4080485/ERR4080485_1.fastq.gz",
             "fastq_md5": "4de269d2b5831e1c5175586af694d21e",
             "host_tax_id": "9606",
             "collection_date": "2019-12-31 12:12:12"}
        ])
        ena_accessor.access()
        self.assertEqual(ena_accessor.included, 0)
        self.assertEqual(ena_accessor.excluded, 3)
