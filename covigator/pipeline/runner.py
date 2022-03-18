from logzero import logger
import time
import subprocess
from covigator.exceptions import CovigatorPipelineError


def run_command(command, temporary_folder):
    start = time.time()
    p = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=temporary_folder, shell=True)
    try:
        stdoutdata, stderrdata = p.communicate(timeout=1800)    # hard coded timeout to 30 minutes
    except subprocess.TimeoutExpired:
        raise CovigatorPipelineError("Timeout in pipeline: {}".format(command))
    logger.info("Finished in {} secs command: '{}'".format(time.time() - start, command))
    if p.returncode != 0:
        error_message = decode(stderrdata)
        logger.error(error_message)
        raise CovigatorPipelineError("Error executing pipeline command: {}\n{}".format(command, error_message))


def decode(data):
    return data.decode("utf8")
