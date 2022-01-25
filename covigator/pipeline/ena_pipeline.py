import os
from pathlib import Path
from logzero import logger
from covigator.configuration import Configuration
from covigator.pipeline.runner import run_command


class Pipeline:

    def __init__(self, config: Configuration):
        self.config = config
        assert self.config.reference_genome is not None and os.path.exists(self.config.reference_genome), \
            "Please configure the reference genome in the variable {}".format(self.config.ENV_COVIGATOR_REF_FASTA)

    def run(self, run_accession: str, fastq1: str, fastq2: str = None):

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        sample_data_folder = Path(fastq1).parent

        logger.info("Sample data folder: {}".format(sample_data_folder))
        # TODO: this can be simplified to the sample_data_folder
        output_vcf = os.path.join(
            self.config.storage_folder, run_accession,
            "{name}.lofreq.normalized.annotated.vcf.gz".format(name=run_accession))
        output_qc = os.path.join(
            self.config.storage_folder, run_accession, "{name}.fastp_stats.json".format(name=run_accession))
        output_vertical_coverage = os.path.join(
            self.config.storage_folder, run_accession, "{name}.depth.tsv".format(name=run_accession))
        output_horizontal_coverage = os.path.join(
            self.config.storage_folder, run_accession, "{name}.coverage.tsv".format(name=run_accession))

        if not os.path.exists(output_vcf) \
                or not os.path.exists(output_qc) \
                or not os.path.exists(output_horizontal_coverage) \
                or not os.path.exists(output_vertical_coverage) \
                or self.config.force_pipeline:

            command = "{nextflow} run {workflow} " \
                      "--fastq1 {fastq1} {fastq2} --output {output_folder} --name {name} " \
                      "--cpus {cpus} --memory {memory} " \
                      "-profile conda -offline -work-dir {work_folder} -with-trace {trace_file}".format(
                nextflow=self.config.nextflow,
                fastq1=fastq1,
                fastq2="--fastq2 " + fastq2 if fastq2 else "",
                output_folder=sample_data_folder,
                name=run_accession,
                work_folder=self.config.temp_folder,
                workflow=self.config.workflow,
                trace_file=os.path.join(sample_data_folder, "nextflow_traces.txt"),
                cpus=self.config.workflow_cpus,
                memory=self.config.workflow_memory
            )
            run_command(command, sample_data_folder)

        return output_vcf, output_qc, output_vertical_coverage, output_horizontal_coverage
