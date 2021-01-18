#!/usr/bin/env python

from argparse import ArgumentParser
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import List

from logzero import logger

import config as cfg

class Pipeline:

    def __init__(self, fastqs : List[str]):
        self.fastqs = fastqs

    def run(self):
        logger.info("Processing {}".format(self.fastqs))
        time.sleep(10)

        cmds = cfg.cmds
        refs = cfg.references
        fq1 = self.fastqs[0]
        fq2 = self.fastqs[1] if len(self.fastqs) > 1 else None
        fq_path = Path(self.fastqs[0]).parent
        virus_db_path = refs["virus_db"]
        bwa_ref = refs["novo_fasta"]
        
        # Creating a temporary folder for the intermediate results
        # and copy the final VCF files to FASTQ folder
        with tempfile.TemporaryDirectory() as tmpdir:
            sam_file = os.path.join(tmpdir, "Aligned.out.sam")
            bam_file = os.path.join(tmpdir, "Aligned.out.bam")
            vcf_file = os.path.join(fq_path, "pileup.vcf")
            snpeff_vcf_file = os.path.join(fq_path, "snpeff.vcf")

            if fq2:
                cmd_align = "{} mem {} {} {} > {}".format(cmds["bwa"], bwa_ref, fq1, fq2, sam_file)
            else:
                cmd_align = "{} mem {} {} > {}".format(cmds["bwa"], bwa_ref, fq1, sam_file)
            cmd_sambam = "{} sort -O bam -o {} {}".format(cmds["samtools"], bam_file, sam_file)

            # Currently the pipeline supports both mpileup and lofreq
            cmd_pileup = "{0} mpileup -E -d 0 -A -f {1} {2} | {0} call -mv --ploidy 1 -Ov -o {3}".format(cmds["bcftools"], refs["novo_fasta"], bam_file, vcf_file)
#            cmd_lofreq = "{0} call -s 0.01 -q 0 -Q 0 --no-default-filter -f {1} -o {2} {3}".format(cmds["lofreq"], refs["novo_fasta"], vcf_file, bam_file)
            cmd_snpeff = "{} ann -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -fi {} SARS-COV2 {} > {}".format(cmds["snpeff"], refs["novo_bed"], vcf_file, snpeff_vcf_file)

            for cmd in [
                    cmd_align,
                    cmd_sambam,
 #                   cmd_lofreq,
                    cmd_pileup,
                    cmd_snpeff
            ]:
                logger.info("Executing CMD: {}".format(cmd))
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmpdir, shell=True)
                (stdoutdata, stderrdata) = p.communicate()
                r = p.returncode
                if r != 0:
                    logger.error(stderrdata)
                    sys.exit(1)

        return snpeff_vcf_file

if __name__ == "__main__":
    # Implementation of Argument Parser for easy testing of the pipeline
    # Maybe implement Unit Test in future
    parser = ArgumentParser(description="Run Pipeline for testing")
    parser.add_argument("-i", "--input", dest="input_fastqs", nargs="+", help="Specify space separated list of paired-end fastqs to process.", required=True)

    args = parser.parse_args()
    pipeline = Pipeline(args.input_fastqs)
    vcf_file = pipeline.run()
