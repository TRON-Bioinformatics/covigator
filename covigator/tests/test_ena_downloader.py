from unittest import TestCase
from covigator.download.ena_downloader import EnaDownloader
from dask.distributed import Client


class EnaDownloaderTests(TestCase):

    def test_download(self):
        sarscov2_tax_id = "2697049"
        homo_sapiens_tax_id = "9606"
        dask_client = Client(n_workers=2, threads_per_worker=1)
        downloader = EnaDownloader(
            tax_id=sarscov2_tax_id, host_tax_id=homo_sapiens_tax_id, dask_client=dask_client, data_folder=None)
        downloader.download()
