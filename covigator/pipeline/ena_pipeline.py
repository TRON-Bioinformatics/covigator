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
        output_vcf = os.path.join(
            self.config.storage_folder, run_accession,
            "{name}.lofreq.normalized.annotated.vcf.gz".format(name=run_accession))
        output_qc = os.path.join(
            self.config.storage_folder, run_accession,
            "{name}.fastp_stats.json".format(name=run_accession))

        if not os.path.exists(output_vcf) or not os.path.exists(output_qc) or self.config.force_pipeline:

            command = "{nextflow} run {workflow} " \
                      "{tronflow_bwa} " \
                      "{tronflow_bam_preprocessing} "\
                      "{tronflow_variant_normalization} " \
                      "--fastq1 {fastq1} {fastq2} --output {output_folder} --name {name} " \
                      "--cpus {cpus} --memory {memory}" \
                      "-profile conda -offline -work-dir {work_folder} -with-trace {trace_file}".format(
                nextflow=self.config.nextflow,
                fastq1=fastq1,
                fastq2="--fastq2 " + fastq2 if fastq2 else "",
                output_folder=self.config.storage_folder,
                name=run_accession,
                work_folder=self.config.temp_folder,
                workflow=self.config.workflow,
                tronflow_bwa="--tronflow_bwa {}".format(self.config.tronflow_bwa) if self.config.tronflow_bwa else "",
                tronflow_bam_preprocessing="--tronflow_bam_preprocessing {}".format(
                    self.config.tronflow_bam_preprocessing) if self.config.tronflow_bam_preprocessing else "",
                tronflow_variant_normalization="--tronflow_variant_normalization {}".format(
                    self.config.tronflow_variant_normalization) if self.config.tronflow_variant_normalization else "",
                trace_file=os.path.join(sample_data_folder, "nextflow_traces.txt"),
                cpus=self.config.workflow_cpus,
                memory=self.config.workflow_memory
            )
            self._run_command(command, sample_data_folder)

        return output_vcf, output_qc

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
