import os
import time
from pathlib import Path
import subprocess
from logzero import logger
from covigator.configuration import Configuration
from covigator.exceptions import CovigatorPipelineError


class Pipeline:

    def __init__(self, config: Configuration):
        self.config = config

    def run(self, run_accession: str, fastq1: str, fastq2: str = None):

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        sample_data_folder = Path(fastq1).parent

        logger.info("Sample data folder: {}".format(sample_data_folder))
        snpeff_vcf_file_gz = os.path.join(
            self.config.storage_folder, run_accession,
            "{name}.lofreq.normalized.annotated.vcf.gz".format(name=run_accession))

        command = "{nextflow} run {workflow} -r {workflow_version} " \
                  "--fastq1 {fastq1} {fastq2} --output {output_folder} --name {name} " \
                  "-profile conda -offline -work-dir {work_folder}".format(
            nextflow=self.config.nextflow,
            fastq1=fastq1,
            fastq2="--fastq2 " + fastq2 if fastq2 else "",
            output_folder=self.config.storage_folder,
            name=run_accession,
            work_folder=self.config.temp_folder,
            workflow_version=self.config.workflow_version,
            workflow=self.config.workflow)

        self._run_command(command, sample_data_folder)

        return snpeff_vcf_file_gz

    def _run_command(self, command, temporary_folder):
            start = time.time()
            p = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=temporary_folder, shell=True)
            stdoutdata, stderrdata = p.communicate()
            logger.info("Finished in {} secs command: '{}'".format(time.time() - start, command))
            if p.returncode != 0:
                error_message = self._decode(stderrdata)
                logger.error(error_message)
                raise CovigatorPipelineError("Error executing pipeline command: {}\n{}".format(command, error_message))

    def _decode(self, data):
        return data.decode("utf8")
