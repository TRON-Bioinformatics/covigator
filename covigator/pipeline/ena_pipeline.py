import os
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
        fq_path = Path(fastq1).parent
        bwa_ref = self.config.reference_genome
        
        # Creating a temporary folder for the intermediate results
        # and copy the final VCF files to FASTQ folder
        with tempfile.TemporaryDirectory() as tmpdir:
            logger.info("Temporary folder: {}".format(tmpdir))
            bam_file = os.path.join(tmpdir, "aligned.out.bam")
            vcf_file = os.path.join(tmpdir, "pileup.vcf")
            snpeff_vcf_file_gz = os.path.join(fq_path, "snpeff.vcf.gz")

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

            self._run_commands([cmd_align, cmd_pileup, cmd_snpeff], tmpdir)

        return snpeff_vcf_file_gz

    def _run_commands(self, commands, temporary_folder):
        for command in commands:
            logger.info("Executing: {}".format(command))
            p = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=temporary_folder, shell=True)
            stdoutdata, stderrdata = p.communicate()
            if p.returncode != 0:
                logger.error(stderrdata)
                raise CovigatorPipelineError("Error executing pipeline command: {}".format(command))
