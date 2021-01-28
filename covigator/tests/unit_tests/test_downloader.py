from unittest import TestCase
from covigator.model import EnaRun
from covigator.processor.downloader import Downloader, CovigatorMD5CheckSumError
import os


class DownloaderTest(TestCase):

    def setUp(self) -> None:
        self.downloader = Downloader()

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
            fastq_ftp="https://tron-mainz.de/wp-content/uploads/2020/07/TRON_Logo_Science.svg",
            fastq_md5="ca72bacad0dfcf665df49bfc53cc8b60"
        )
        self.downloader.download(ena_run)
        run_storage_folder = os.path.join(self.downloader.storage_folder, run_accession)
        self.assertTrue(os.path.exists(run_storage_folder))
        self.assertTrue(os.path.exists(os.path.join(run_storage_folder, "TRON_Logo_Science.svg")))

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
