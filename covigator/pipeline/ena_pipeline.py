import os
from dataclasses import dataclass
from pathlib import Path
from logzero import logger
from covigator.configuration import Configuration
from covigator.pipeline.runner import run_command


@dataclass
class PipelineResult:
    # VCF
    lofreq_vcf: str
    ivar_vcf:  str
    gatk_vcf: str
    bcftools_vcf: str
    # QC
    fastp_qc: str
    vertical_coverage: str
    horizontal_coverage: str
    deduplication_metrics: str
    # pangolin results
    lofreq_pangolin: str
    ivar_pangolin: str
    gatk_pangolin: str
    bcftools_pangolin: str


class Pipeline:

    def __init__(self, config: Configuration):
        self.config = config
        assert self.config.reference_genome is not None and os.path.exists(self.config.reference_genome), \
            "Please configure the reference genome in the variable {}".format(self.config.ENV_COVIGATOR_REF_FASTA)

    def run(self, run_accession: str, fastq1: str, fastq2: str = None) -> PipelineResult:

        logger.info("Processing {} and {}".format(fastq1, fastq2))
        sample_data_folder = Path(fastq1).parent

        logger.info("Sample data folder: {}".format(sample_data_folder))

        lofreq_vcf = os.path.join(sample_data_folder, "{name}.lofreq.vcf.gz".format(name=run_accession))
        ivar_vcf = os.path.join(sample_data_folder, "{name}.ivar.vcf.gz".format(name=run_accession))
        gatk_vcf = os.path.join(sample_data_folder, "{name}.gatk.vcf.gz".format(name=run_accession))
        bcftools_vcf = os.path.join(sample_data_folder, "{name}.bcftools.vcf.gz".format(name=run_accession))

        lofreq_pangolin = os.path.join(sample_data_folder, "{name}.lofreq.pangolin.csv".format(name=run_accession))
        ivar_pangolin = os.path.join(sample_data_folder, "{name}.ivar.pangolin.csv".format(name=run_accession))
        gatk_pangolin = os.path.join(sample_data_folder, "{name}.gatk.pangolin.csv".format(name=run_accession))
        bcftools_pangolin = os.path.join(sample_data_folder, "{name}.bcftools.pangolin.csv".format(name=run_accession))

        output_qc = os.path.join(sample_data_folder, "{name}.fastp_stats.json".format(name=run_accession))
        output_vertical_coverage = os.path.join(sample_data_folder, "{name}.depth.tsv".format(name=run_accession))
        output_horizontal_coverage = os.path.join(sample_data_folder, "{name}.coverage.tsv".format(name=run_accession))
        deduplication_metrics = os.path.join(
            sample_data_folder, "{name}.deduplication_metrics.txt".format(name=run_accession))

        if not os.path.exists(lofreq_vcf) \
                or not os.path.exists(output_qc) \
                or not os.path.exists(output_horizontal_coverage) \
                or not os.path.exists(output_vertical_coverage) \
                or self.config.force_pipeline:

            command = "{nextflow} run {workflow} " \
                      "--fastq1 {fastq1} {fastq2} " \
                      "--output {output_folder} " \
                      "--name {name} " \
                      "--low_frequency_variant_threshold {af_low_frequency_thr} " \
                      "--subclonal_variant_threshold {af_subclonal_thr} " \
                      "--cpus {cpus} " \
                      "--memory {memory} " \
                      "--skip_bcftools " \
                      "--skip_gatk " \
                      "-profile conda " \
                      "-offline " \
                      "-work-dir {work_folder} " \
                      "-with-trace {trace_file}".format(
                nextflow=self.config.nextflow,
                fastq1=fastq1,
                fastq2="--fastq2 " + fastq2 if fastq2 else "",
                af_low_frequency_thr=self.config.low_coverage_threshold,
                af_subclonal_thr=self.config.subclonal_threshold,
                output_folder=sample_data_folder,
                name=run_accession,
                work_folder=self.config.temp_folder,
                workflow=self.config.workflow,
                trace_file=os.path.join(sample_data_folder, "nextflow_traces.txt"),
                cpus=self.config.workflow_cpus,
                memory=self.config.workflow_memory)
            run_command(command, sample_data_folder)

        return PipelineResult(
            # VCF
            lofreq_vcf=lofreq_vcf,
            ivar_vcf=ivar_vcf,
            gatk_vcf=gatk_vcf,
            bcftools_vcf=bcftools_vcf,
            # QC
            fastp_qc=output_qc,
            vertical_coverage=output_vertical_coverage,
            horizontal_coverage=output_horizontal_coverage,
            deduplication_metrics=deduplication_metrics,
            # pangolin
            lofreq_pangolin=lofreq_pangolin,
            ivar_pangolin=ivar_pangolin,
            gatk_pangolin=gatk_pangolin,
            bcftools_pangolin=bcftools_pangolin
        )
