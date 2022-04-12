import os
import pathlib
import urllib.request as request
import hashlib
from logzero import logger
from covigator.configuration import Configuration
from covigator.database.model import SampleEna, SEPARATOR
import re

FTP_PROTOCOL = "ftp://"
PROTOCOL_REGEX = re.compile('[a-zA-Z0-9]+://')


class CovigatorMD5CheckSumError(Exception):
    pass


class Downloader:

    def __init__(self, config: Configuration):
        self.storage_folder = config.storage_folder
        pathlib.Path(self.storage_folder).mkdir(parents=True, exist_ok=True)
        assert os.path.exists(self.storage_folder), "Storage folder does not exist"

    def download(self, sample_ena: SampleEna) -> str:
        assert sample_ena.run_accession is not None, "A run accession is required as it is part of the file path"
        assert sample_ena.fastq_ftp is not None, "Cannot download empty URLs"
        assert sample_ena.fastq_md5 is not None, "Cannot do MD5 checks without MD5"
        fastqs = sample_ena.get_fastqs_ftp()
        md5s = sample_ena.get_fastqs_md5()
        paths = []
        for url, md5 in zip(fastqs, md5s):
            local_path = self._download_url(sample_ena, url, md5)
            paths.append(local_path)
        logger.info("Downloaded {}: {}".format(sample_ena.run_accession, sample_ena.fastq_ftp))
        return SEPARATOR.join(paths)

    def _download_url(self, sample: SampleEna, url, md5):
        """
        This method streams a file (probably large) from a URL and stores it in the local file system
        """
        if not re.match(PROTOCOL_REGEX, url):
            url = FTP_PROTOCOL + url
        local_filename = url.split('/')[-1]
        # NOTE: sample folder date/run_accession
        local_folder = sample.get_sample_folder(self.storage_folder)
        local_full_path = os.path.join(local_folder, local_filename)
        # avoids downloading the same files over and over
        if not os.path.exists(local_full_path):
            pathlib.Path(local_folder).mkdir(parents=True, exist_ok=True)
            request.urlretrieve(url, local_full_path)
            self._check_md5(local_full_path, md5)

        return local_full_path

    def _check_md5(self, filepath, md5):
        """
        This method reads a file (probably large) in chunks and computes the MD5 hash iteratively
        """
        file_hash = hashlib.md5()
        with open(filepath, "rb") as f:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
                file_hash.update(chunk)
        if file_hash.hexdigest() != md5:
            raise CovigatorMD5CheckSumError("Failed MD5 check for file {} and MD5 {}".format(filepath, md5))
