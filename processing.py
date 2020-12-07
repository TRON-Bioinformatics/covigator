#!/usr/bin/env python

# Scripts
# Minion scripts for corona virus pipeline
# /projects/SARS-CoV-2/minion_gridion.sh

# RNA-Seq script for corona virus pipeline
# /projects/SARS-CoV-2/documentation_THBU.sh


# Data
# RNA-Seq data for corona virus pipeline
# /scratch/info/projects/SARS-CoV-2/Novoalign/fastqs/

# Minion data for corona virus pipeline
# /scratch/info/projects/SARS-CoV-2/sra2

# GISAID assembly data for Coronavirus
# /scratch/info/projects/SARS-CoV-2/gisaid/


from argparse import ArgumentParser
import os
from shutil import copy


import config as cfg

import misc.io_methods as IOMethods
from misc.queue import Queue

class Processing(object):
    def __init__(self, input_paths, working_dir):
        self.working_dir = os.path.abspath(working_dir)
        #self.logger = Logger(os.path.join(self.working_dir, "easyfuse_processing.log"))
        IOMethods.create_folder(self.working_dir)
        copy(os.path.join(cfg.module_dir, "config.py"), working_dir)

        self.input_paths = [os.path.abspath(file) for file in input_paths]

    def run(self):

        fastqs = IOMethods.get_fastq_files(self.input_paths)
        left, right, sample_id = IOMethods.pair_fastq_files(fastqs)

        for i, _ in enumerate(left):
            if len(left) == len(right):
                self.logger.info("Processing Sample ID: {} (paired end)".format(sample_id[i]))
                self.logger.info("Sample 1: {}".format(left[i]))
                self.logger.info("Sample 2: {}".format(right[i]))
                self.execute_pipeline(left[i], right[i], sample_id[i])

    def execute_pipeline(self, fq1, fq2, sample_id):

        cmds = cfg.cmds

        refs = cfg.references

        # build virus reference DB
        virus_db_path = refs["virus_db"]
        bwa_ref = refs["novo_fasta"]

        #IOMethods.create_folder(virus_db_path)

        #download_virus_db(virus_db_path)

        # download GISAID

        # download SRA

        # download reference nt

        # do pairwise alignment

        # download NGS reference

        # do BWA alignment

        alignment_path = os.path.join(self.working_dir, "alignment")
        sam_file = os.path.join(alignment_path, "Aligned.out.sam")
        bam_file = os.path.join(alignment_path, "Aligned.out.bam")
        
        cmd_align = "{} mem {} {} {} > {}".format(cmds["bwa"], bwa_ref, fq1, fq2, sam_file)
        cmd_sambam = "{} sort -O bam -o {} {}".format(cmds["samtools"], bam_file, sam_file)
        
        # do mpileup

        cmd_pileup = "{0} mpileup -E -d 0 -A -f {1} {2} | {0} call -mv --ploidy 1 -Ov -o {3}".format(cmds["bcftools"], refs["novo_fasta"], bam_file, vcf_file)

        # do snpeff

        cmd_snpeff = "{} ann -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -fi {} SARS-COV2 {} > {}".format(cmds["snpeff"], refs["novo_bed"], vcf_file, snpeff_vcf_file)

        exe_tools = [
            "BWA",
            "PileupVCF",
            "SNPeff",
            "Collect"
        ]

        exe_cmds = [
            " && ".join([cmd_align, cmd_sambam]),
            " && ".join([cmd_pileup, cmd_filter]),
            cmd_snpeff,
            cmd_gather
        ]

        exe_path = [
            alignment_path,
            res_path,
            res_path,
            res_path
        ]

        for i, tool in enumerate(exe_tools, 0):
            if tool in tools:
                dependency = []
                if tool in state_tools:
                    self.logger.info("Skipping {0} as it looks like a previous run finished successfully. Results should be in {1}".format(tool, exe_path[i]))
                    continue

                self.logger.info("Submitting {} run to slurm".format(tool))
                cpu = str(cfg.resources[tool.lower()]["cpu"])
                mem = cfg.resources[tool.lower()]["mem"]

                exe_cmds[i] = exe_cmds[i].replace("waiting_for_cpu_number", cpu)

                dependency = Queueing.get_jobs_by_name(sample_id, cfg.pipeline_name)
                uid = "-".join([cfg.pipeline_name, tool, sample_id])
                cmd = " && ".join([exe_cmds[i], cmd_samples + tool])
                self.submit_job(uid, cmd, cpu, mem, exe_path[i], dependency, "")


    def submit_job(self, uid, cmd, cores, mem_usage, job_dir, dependencies, mail):
        """This function submits a job with the corresponding resources to slurm."""
        already_running = Queueing.get_jobs_by_name(uid)
        if not already_running:
            module_file = os.path.join(cfg.module_dir, "build_env.sh")
            que_sys = cfg.queueing_system
            for i, cmd_split in enumerate(cmd.split(" && ")):
                if not que_sys in ["slurm", "pbs"]:
                    cmd_split = cmd_split.split(" ")
                dependencies.extend(Queueing.get_jobs_by_name("{0}_CMD{1}".format(uid, i - 1)))
                Queueing.submit(
                    "{0}_CMD{1}".format(uid, i), 
                    cmd_split, 
                    cores, 
                    mem_usage, 
                    job_dir, 
                    dependencies, 
                    cfg.partition, 
                    cfg.user, 
                    cfg.time_limit, 
                    mail, 
                    module_file, 
                    que_sys
                )
                time.sleep(0.5)

        else:
            self.logger.info(uid + " already running!")
            

def download_virus_db(db_path):
    cmds = []
    download_cmd_1 = "wget -N https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
    download_cmd_2 = "wget -N https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"

    # inflate and write to 1 file
    extract_cmd = "zcat viral.*.genomic.fna.gz > aux.viruses.fna"

    # changing 'chromosome' names from 'gi|167006425|ref|NC_010314.1|' to NC_010314.1
    sed_cmd = "sed 's/^>gi|[0-9]*|ref|/>/g;s/[;'\'']/_/g;s/|//g' aux.viruses.fna > all.viruses.fna"

    ids_cmd = "grep -oP '[A-Z]{2,2}_[0-9]*\.[0-9]{1,1}' all.viruses.fna > all.viruses.IDs.txt"

    cmds.append(("Download 1", download_cmd_1))
    cmds.append(("Download 2", download_cmd_2))
    cmds.append(("Zcat", extract_cmd))
    cmds.append(("Replace IDs", sed_cmd))
    cmds.append(("Extract IDs", ids_cmd))


    q = Queue()
    for name, cmd in cmds:
        print("{}: {}".format(name, cmd))
        q.submit_nonqueue(cmd, db_path)

def main():
    parser = ArgumentParser(description="Automatically runs analysis pipeline")
    parser.add_argument('-i', '--input', dest='input_paths', nargs='+', help='Specify full path of the fastq folder to process.', required=True)
    parser.add_argument("-w", "--working_dir", dest="working_dir", help="Specify working directory")

    args = parser.parse_args()
    
    proc = Processing(args.input_paths, args.working_dir)

    proc.run()
    

if __name__ == "__main__":
    main()
