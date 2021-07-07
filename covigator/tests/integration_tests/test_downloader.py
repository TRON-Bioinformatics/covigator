from unittest import TestCase

from covigator.configuration import Configuration
from covigator.database.model import SampleEna
from covigator.pipeline.downloader import Downloader
import os


class DownloaderTest(TestCase):

    def setUp(self) -> None:
        self.downloader = Downloader(config=Configuration())

    def test_download_ena_run(self):
        run_accession = "TEST12345"
        ena_run = SampleEna(
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
