import os
from pathlib import Path
import subprocess
import tempfile
from logzero import logger


class CovigatorPipelineError(Exception):
    pass


class Pipeline:

    # TODO: make this configurable through environment variables
    commands = {
        "samtools": "/code/samtools/1.9/samtools",
        "bwa": "/code/bwa/0.7.17/bwa",
        "bcftools": "/code/bcftools/1.9/bcftools",
        "snpeff": "java -jar /code/snpEFF_latest_core/snpEff/snpEff.jar",
    }

    # TODO: automate the download of references and configure location through environment variable
    references = {
        "novo_fasta": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3.fa",
        "novo_bed": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3_gp02_Sgene.bed"
    }

    def run(self, fastq1: str, fastq2: str = None):

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        fq_path = Path(fastq1).parent
        bwa_ref = self.references["novo_fasta"]
        
        # Creating a temporary folder for the intermediate results
        # and copy the final VCF files to FASTQ folder
        with tempfile.TemporaryDirectory() as tmpdir:
            logger.info("Temporary folder: {}".format(tmpdir))
            sam_file = os.path.join(tmpdir, "aligned.out.sam")
            bam_file = os.path.join(tmpdir, "aligned.out.bam")
            vcf_file = os.path.join(fq_path, "pileup.vcf")
            snpeff_vcf_file = os.path.join(fq_path, "snpeff.vcf")

            if fastq2:
                cmd_align = "{} mem {} {} {} > {}".format(self.commands["bwa"], bwa_ref, fastq1, fastq2, sam_file)
            else:
                cmd_align = "{} mem {} {} > {}".format(self.commands["bwa"], bwa_ref, fastq1, sam_file)
            cmd_sambam = "{} sort -O bam -o {} {}".format(self.commands["samtools"], bam_file, sam_file)

            # Currently the pipeline supports both mpileup and lofreq
            cmd_pileup = "{0} mpileup -E -d 0 -A -f {1} {2} | {0} call -mv --ploidy 1 -Ov -o {3}".format(
                self.commands["bcftools"], self.references["novo_fasta"], bam_file, vcf_file)
            cmd_snpeff = "{} ann -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -fi {} " \
                         "SARS-COV2 {} > {}".format(
                self.commands["snpeff"], self.references["novo_bed"], vcf_file, snpeff_vcf_file)

            self._run_commands([cmd_align, cmd_pileup, cmd_sambam, cmd_snpeff], tmpdir)

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
