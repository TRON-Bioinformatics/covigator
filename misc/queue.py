#!/usr/bin/env python

import os
import sys
import time
import subprocess
from argparse import ArgumentParser

class Queue(object):
    
    @staticmethod
    def submit_nonqueue(cmd, output_results_folder):
#        print output_results_folder
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_results_folder, shell=True)
        (stdoutdata, stderrdata) = p.communicate()
        r = p.returncode
        if r != 0:
            print(stderrdata)
            sys.exit(1)

    @staticmethod
    def submit_pbs(job_name, cmd, cores, mem_usage, output_results_folder, dependencies):
#        job_name = "-".join([module,sample_id])
        
        p = subprocess.Popen('qsub', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        (output, input) = (p.stdout, p.stdin)

    # Customize your options here
#        walltime = "999:00:00"
        processors = "nodes=1:ppn=" + str(cores) + ",mem=" + str(mem_usage) + "gb,vmem=" + str(mem_usage) + "gb"
        error_file = os.path.join(output_results_folder, "error.log")
        output_file = os.path.join(output_results_folder, "output.log")

        job_string = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=999:00:00
#PBS -l %s
#PBS -W depend=afterok:%s
#PBS -d %s
#PBS -e %s
#PBS -o %s

%s""" % (job_name, walltime, processors, ":".join(dependencies), output_results_folder, error_file, output_file, cmd)

        # Send job_string to qsub
        input.write(job_string)
        input.close()
#        log(output_results_folder,job_string,"STATUS")
        # Print your job and the response to the screen
        print(job_string)
        #    print output.read()
        
        time.sleep(1)


    @staticmethod
    def submit_slurm(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, partitions, reserved=""):
        '''This function submits a predefined job with specific SBATCH parameters to the Slurm workload manager system.'''
#        job_name = "-".join([module,sample_id])

        p = subprocess.Popen('sbatch', stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        (output, input) = (p.stdout, p.stdin)

        # Customize your options here
        partitions_str = ""
        if reserved:
            partitions_str = "%s -x %s" % (partitions, reserved)
        else:
            partitions_str = "%s" % (partitions)
        depend = ""
        if len(dependencies) != 0:
            depend = "#SBATCH --dependency=afterok:%s" % ":".join(dependencies)
        else:
            depend = ""
        error_file = os.path.join(output_results_folder, "error.log")
        output_file = os.path.join(output_results_folder, "output.log")

        job_string = """#!/bin/bash
#SBATCH -J %s
#SBATCH -p %s
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --cpus-per-task=%s
#SBATCH --mem=%s
#SBATCH --time=30-00:00:00
%s
#SBATCH --workdir=%s
#SBATCH --error=%s
#SBATCH --output=%s

. /etc/profile.d/modules.sh
module load software/Java/openJDK-1.8
module load software/python/python-2.7.9

srun %s""" % (job_name, partitions_str, cores, int(mem_usage)*1000, depend, output_results_folder, error_file, output_file, cmd)

        # Send job_string to qsub
        input.write(job_string)
        input.close()
#        log(output_results_folder,job_string,"STATUS")
        # Print your job and the response to the screen
        outf = open(os.path.join(output_results_folder, job_name + ".sbatch"), "w")
        print(job_string)
        outf.write(job_string)
        outf.close()
        #    print output.read()

#        time.sleep(1800)
#. /etc/profile.d/modules.sh
#module load software/Java/openJDK-1.8

        

def main():
    parser = ArgumentParser(description='Handle data streams')
    parser.add_argument('-j', '--job-name', dest='job_name', help='Specify the job name.', required=True )
    parser.add_argument('-s', '--script', dest='script', help='Specify the script to process.', required=True )
    parser.add_argument('-c', '--cores', dest='cores', type=int, help='Specify the number of cores to run the script with.', required=True )
    parser.add_argument('-m', '--memory', dest='memory', type=int, help='Specify the amount of memory to run your script with.', required=True)
    parser.add_argument('-o', '--output', dest='output', help='Select the output folder.', required=True)
    parser.add_argument('-p', '--partitions', dest='partitions', help='Select the slurm partitions for this job.', default='allNodes')
    args = parser.parse_args()

    job_name = args.job_name # "STAR-0116_094_RL"
    cmd = args.script # "cd /home/sorn/galaxy-central-vm/scripts/api/demultiplexing_bot/scratch/test123 && /kitty/code/STAR_2.4.2a/bin/Linux_x86_64/STAR --genomeDir /kitty/data/human/STAR_idx/2.3.0/human/genome_SJ/ --genomeLoad NoSharedMemory --readFilesCommand 'gzip -d -c -f' --readFilesIn /kitty/data/seq/140919_SN138_0289_BC4DJNACXX/Project_921964988/Sample_0116_094_RL/0116_094_RL_GTGAAA_L005_R1_001.fastq.gz /kitty/data/seq/140919_SN138_0289_BC4DJNACXX/Project_921964988/Sample_0116_094_RL/0116_094_RL_GTGAAA_L005_R2_001.fastq.gz --outSAMmode Full --outSAMattributes Standard --outSAMunmapped None --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.02 --runThreadN 6 --outSAMtype BAM SortedByCoordinate && touch aligned.pid"
    cores = args.cores # 6
    mem_usage = args.memory # 40
    output_results_folder = args.output # "/home/sorn/galaxy-central-vm/scripts/api/demultiplexing_bot/scratch/test123"
    dependencies = ""
    Queue.submit_slurm(job_name, cmd, cores, mem_usage, output_results_folder, dependencies, args.partitions)

if __name__ == '__main__':
    main()
