import os
import time
from pathlib import Path
import subprocess
import tempfile
from logzero import logger
from covigator.configuration import Configuration
from covigator.exceptions import CovigatorPipelineError


class Pipeline:

    def __init__(self, config: Configuration):
        self.config = config

    def run(self, fastq1: str, fastq2: str = None):

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        sample_data_folder = Path(fastq1).parent
        bwa_ref = self.config.reference_genome

        logger.info("Sample data folder: {}".format(sample_data_folder))
        bam_file = os.path.join(sample_data_folder, "aligned.out.bam")
        vcf_file = os.path.join(sample_data_folder, "pileup.vcf")
        snpeff_vcf_file_gz = os.path.join(sample_data_folder, "snpeff.vcf.gz")

        if fastq2:
            cmd_align = "{} mem {} {} {} | {} sort -o {} -".format(
                self.config.bwa, bwa_ref, fastq1, fastq2, self.config.samtools, bam_file)
        else:
            cmd_align = "{} mem {} {} | {} sort -o {} -".format(
                self.config.bwa, bwa_ref, fastq1, self.config.samtools, bam_file)

        cmd_pileup = "{0} mpileup -E -d 0 -A -f {1} {2} | {0} call -mv --ploidy 1 -Ov -o {3}".format(
            self.config.bcftools, bwa_ref, bam_file, vcf_file)

        cmd_snpeff = "{java} -jar {snpeff} ann " \
                     "-noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein " \
                     "SARS-COV2 {input_vcf} | {bgzip} -c > {annotated_vcf} && {tabix} -p vcf {annotated_vcf}"\
            .format(java=self.config.java, snpeff=self.config.snpeff, input_vcf=vcf_file, bgzip=self.config.bgzip,
                    annotated_vcf=snpeff_vcf_file_gz, tabix=self.config.tabix)

        self._run_commands([cmd_align, cmd_pileup, cmd_snpeff], sample_data_folder)
        self.delete_files([bam_file, vcf_file])

        return snpeff_vcf_file_gz

    def _run_commands(self, commands, temporary_folder):
        for command in commands:
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

    def delete_files(self, files):
        for f in files:
            os.remove(f)
