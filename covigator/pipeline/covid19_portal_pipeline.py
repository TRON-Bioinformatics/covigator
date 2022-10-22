#!/usr/bin/env python
import os
from dataclasses import dataclass
from covigator.configuration import Configuration
from covigator.database.model import SampleCovid19Portal
from covigator.pipeline.runner import run_command


@dataclass
class Covid19PortalPipelineResult:
    vcf_path: str
    fasta_path: str
    pangolin_path: str


class Covid19PortalPipeline:

    def __init__(self, config: Configuration):
        self.config = config

    def run(self, sample: SampleCovid19Portal) -> Covid19PortalPipelineResult:
        # NOTE: sample folder date/run_accession
        sample_data_folder = sample.get_sample_folder(self.config.storage_folder)
        output_vcf = os.path.join(sample_data_folder, "{name}.assembly.vcf.gz".format(name=sample.run_accession))
        output_pangolin = os.path.join(sample_data_folder,
                                       "{name}.assembly.pangolin.csv".format(name=sample.run_accession))
        input_fasta = sample.fasta_path

        if not os.path.exists(output_vcf) or self.config.force_pipeline:

            command = "{nextflow} run {workflow} " \
                      "--fasta {fasta} " \
                      "--output {output_folder} " \
                      "--name {name} " \
                      "--cpus {cpus} " \
                      "--memory {memory} " \
                      "-profile conda " \
                      "-offline " \
                      "-work-dir {work_folder} " \
                      "-with-trace {trace_file}".format(
                nextflow=self.config.nextflow,
                fasta=input_fasta,
                output_folder=sample_data_folder,
                name=sample.run_accession,
                work_folder=self.config.temp_folder,
                workflow=self.config.workflow,
                trace_file=os.path.join(sample_data_folder, "nextflow_traces.txt"),
                cpus=self.config.workflow_cpus,
                memory=self.config.workflow_memory
            )
            run_command(command, sample_data_folder)

        return Covid19PortalPipelineResult(
            vcf_path=output_vcf,
            fasta_path=input_fasta,
            pangolin_path=output_pangolin
        )