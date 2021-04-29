import os
import logzero
import logging


VERSION = "0.1.29"


def initialise_logs(logfile):
    if logfile is not None:
        logzero.logfile(logfile, maxBytes=1e6, backupCount=3)
    logzero.loglevel(logging.INFO)
