#!/usr/bin/env python
import os
from dataclasses import dataclass

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from logzero import logger
from covigator.configuration import Configuration
from covigator.database.model import SampleGisaid
from covigator.exceptions import CovigatorExcludedAssemblySequence
from covigator.misc.compression import decompress_sequence
from covigator.pipeline.runner import run_command


@dataclass
class GisaidPipelineResult:
    vcf_path: str
    fasta_path: str
    pangolin_path: str


class GisaidPipeline:

    def __init__(self, config: Configuration):
        self.config = config

    def run(self, sample: SampleGisaid) -> GisaidPipelineResult:
        logger.info("Processing {}".format(sample.run_accession))
        sample_name = sample.run_accession.replace("/", "_").replace(" ", "-").replace("'", "-").replace("$", "")
        # NOTE: sample folder date/run_accession
        sample_data_folder = sample.get_sample_folder(self.config.storage_folder)
        output_vcf = os.path.join(sample_data_folder, "{name}.assembly.vcf.gz".format(name=sample_name))
        output_pangolin = os.path.join(sample_data_folder, "{name}.assembly.pangolin.csv".format(name=sample_name))
        input_fasta = os.path.join(sample_data_folder, "{}.fasta".format(sample_name))

        if not os.path.exists(output_vcf) or self.config.force_pipeline:

            # writes the FASTA sequence into a file
            _, sequence = list(sample.sequence.items())[0]
            os.makedirs(sample_data_folder, exist_ok=True)
            decompressed_sequence = decompress_sequence(sequence)

            # excludes too small sequences
            if len(decompressed_sequence) < self.config.min_sequence_size:
                raise CovigatorExcludedAssemblySequence

            with open(input_fasta, "w+") as output:
                record = SeqRecord(
                    seq=Seq(decompressed_sequence, IUPAC.ambiguous_dna),
                    id=sample_name)
                SeqIO.write(record, output, "fasta")

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
                name=sample_name,
                work_folder=self.config.temp_folder,
                workflow=self.config.workflow,
                trace_file=os.path.join(sample_data_folder, "nextflow_traces.txt"),
                cpus=self.config.workflow_cpus,
                memory=self.config.workflow_memory
            )
            run_command(command, sample_data_folder)

        return GisaidPipelineResult(
            vcf_path=output_vcf,
            fasta_path=input_fasta,
            pangolin_path=output_pangolin
        )
