import os
from pathlib import Path
import subprocess
import tempfile
from logzero import logger

from covigator import ENV_COVIGATOR_BIN_SAMTOOLS, ENV_COVIGATOR_BIN_BWA, ENV_COVIGATOR_BIN_BCFTOOLS, \
    ENV_COVIGATOR_BIN_SNPEFF, ENV_COVIGATOR_BIN_BGZIP, ENV_COVIGATOR_BIN_TABIX, \
    ENV_COVIGATOR_REF_FASTA, ENV_COVIGATOR_REF_BED

class CovigatorPipelineError(Exception):
    pass


class Pipeline:

    commands = {
        "samtools": os.getenv(ENV_COVIGATOR_BIN_SAMTOOLS, "/code/samtools/1.9/samtools"),
        "bwa": os.getenv(ENV_COVIGATOR_BIN_BWA, "/code/bwa/0.7.17/bwa"),
        "bcftools": os.getenv(ENV_COVIGATOR_BIN_BCFTOOLS, "/code/bcftools/1.9/bcftools"),
        "snpeff": os.getenv(ENV_COVIGATOR_BIN_SNPEFF, "/code/snpEFF_latest_core/snpEff/snpEff.jar"),
        "bgzip": os.getenv(ENV_COVIGATOR_BIN_BGZIP, "bgzip"),
        "tabix": os.getenv(ENV_COVIGATOR_BIN_TABIX, "tabix")
    }

    # TODO: automate the download of references and configure location through environment variable

    references = {
        "novo_fasta": os.getenv(ENV_COVIGATOR_REF_FASTA,
                                "/scratch/info/projects/SARS-CoV-2/index/MN908947.3.fa"),
        "novo_bed": os.getenv(ENV_COVIGATOR_REF_BED,
                              "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3_gp02_Sgene.bed")
    }

    def run(self, fastq1: str, fastq2: str = None):

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        fq_path = Path(fastq1).parent
        bwa_ref = self.references["novo_fasta"]
        
        # Creating a temporary folder for the intermediate results
        # and copy the final VCF files to FASTQ folder
        with tempfile.TemporaryDirectory() as tmpdir:
            logger.info("Temporary folder: {}".format(tmpdir))
            bam_file = os.path.join(tmpdir, "aligned.out.bam")
            vcf_file = os.path.join(tmpdir, "pileup.vcf")
            snpeff_vcf_file = os.path.join(tmpdir, "snpeff.vcf")
            snpeff_vcf_file_gz = os.path.join(tmpdir, "snpeff.vcf.gz")
            snpeff_vcf_file_gz_tbi = os.path.join(fq_path, "snpeff.vcf.gz.tbi")

            if fastq2:
                cmd_align = "{} mem {} {} {} | {} sort -o {} -".format(
                    self.commands["bwa"], bwa_ref, fastq1, fastq2, self.commands["samtools"], bam_file)
            else:
                cmd_align = "{} mem {} {} | {} sort -o {} -".format(
                    self.commands["bwa"], bwa_ref, fastq1, self.commands["samtools"], bam_file)

            cmd_pileup = "{0} mpileup -E -d 0 -A -f {1} {2} | {0} call -mv --ploidy 1 -Ov -o {3}".format(
                self.commands["bcftools"], self.references["novo_fasta"], bam_file, vcf_file)

            cmd_snpeff = "java -jar {0} ann -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein "\
                         "SARS-COV2 {1} > {2} && {3} -c {2} > {4} && {5} -p vcf {4} > {6}".format(
                             self.commands["snpeff"], vcf_file, snpeff_vcf_file, self.commands["bgzip"],
                             snpeff_vcf_file_gz, self.commands["tabix"], snpeff_vcf_file_gz_tbi)

            self._run_commands([cmd_align, cmd_pileup, cmd_snpeff], tmpdir)

        return snpeff_vcf_file

    def _run_commands(self, commands, temporary_folder):
        for command in commands:
            logger.info("Executing: {}".format(command))
            p = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=temporary_folder, shell=True)
            stdoutdata, stderrdata = p.communicate()
            if p.returncode != 0:
                logger.error(stderrdata)
                raise CovigatorPipelineError("Error executing pipeline command: {}".format(command))
