from unittest import TestCase
from covigator.model import EnaRun
from covigator.processor.downloader import Downloader, CovigatorMD5CheckSumError
import os


class DownloaderTest(TestCase):

    def setUp(self) -> None:
        self.downloader = Downloader()

    def test_download_ena_run(self):
        run_accession = "TEST12345"
        ena_run = EnaRun(
            run_accession=run_accession,
            fastq_ftp="ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722.fastq.gz;"
                       "ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722_1.fastq.gz;"
                       "ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722_2.fastq.gz",
            fastq_md5="f310c1f2952f07d70aee7464a276cba1;"
                      "11a6108b741b48402001e0bbae16cbaa;"
                      "2dc39ea2a383550c15013435b65489f8"
        )
        self.downloader.download(ena_run)
        run_storage_folder = os.path.join(self.downloader.storage_folder, run_accession)
        self.assertTrue(os.path.exists(run_storage_folder))
        self.assertTrue(os.path.exists(os.path.join(run_storage_folder, "ERR4192722.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(run_storage_folder, "ERR4192722_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(run_storage_folder, "ERR4192722_2.fastq.gz")))

    def test_download_ena_run_with_bad_md5(self):
        run_accession = "TEST12345"
        ena_run = EnaRun(
            run_accession=run_accession,
            fastq_ftp="ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722.fastq.gz;"
                       "ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722_1.fastq.gz;"
                       "ftp.sra.ebi.ac.uk/vol1/fastq/ERR419/002/ERR4192722/ERR4192722_2.fastq.gz",
            fastq_md5="blahblahblah;"
                      "11a6108b741b48402001e0bbae16cbaa;"
                      "2dc39ea2a383550c15013435b65489f8"
        )
        with self.assertRaises(CovigatorMD5CheckSumError):
            self.downloader.download(ena_run)

    def test_download_http(self):
        run_accession = "TEST12346"
        ena_run = EnaRun(
            run_accession=run_accession,
            fastq_ftp="https://commons.wikimedia.org/wiki/File:Google-Logo.svg",
            fastq_md5="74dad1edffff1bc7f08a103ba3565f6e"
        )
        self.downloader.download(ena_run)
        run_storage_folder = os.path.join(self.downloader.storage_folder, run_accession)
        self.assertTrue(os.path.exists(run_storage_folder))
        self.assertTrue(os.path.exists(os.path.join(run_storage_folder, "File:Google-Logo.svg")))

    def test_download_http_without_md5(self):
        run_accession = "TEST12347"
        with self.assertRaises(AssertionError):
            self.downloader.download(EnaRun(
                run_accession=run_accession,
                fastq_ftp="https://commons.wikimedia.org/wiki/File:Google-Logo.svg",
                fastq_md5=None
            ))
        with self.assertRaises(CovigatorMD5CheckSumError):
            self.downloader.download(EnaRun(
                run_accession=run_accession,
                fastq_ftp="https://commons.wikimedia.org/wiki/File:Google-Logo.svg",
                fastq_md5=""
            ))
