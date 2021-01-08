import time
from typing import List

from logzero import logger


class Pipeline:

    def __init__(self, fastqs : List[str]):
        self.fastqs = fastqs

    def run(self):
        logger.info("Processing {}".format(self.fastqs))
        time.sleep(10)
        # TODO: integrate pipeline here
