![CoVigator logo](../figures/CoVigator_logo_txt_reg_no_bg.png "CoVigator logo")

-----------------

# Pipeline

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/covigator-ngs-pipeline)
[![DOI](https://zenodo.org/badge/374669617.svg)](https://zenodo.org/badge/latestdoi/374669617)
[![Run tests](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline/actions/workflows/unit_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline/actions/workflows/unit_tests.yml)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-Nextflow-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

The Covigator pipeline processes SARS-CoV-2 FASTQ or FASTA files into annotated and normalized analysis ready VCF files. 
The pipeline is implemented in the Nextflow framework (Di Tommaso, 2017), it is a stand-alone pipeline that can be
used independently of the CoVigator dashboard and knowledge base. 
The code is open sourced in a repository of its own [https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline).

Although it is configured by default for SARS-CoV-2 it can be employed for the analysis of other microbial organisms 
if the required references are provided.

The result of the pipeline is one or more annotated VCFs with the list of SNVs and indels ready for analysis.

**Table of Contents**

1. [Two pipelines in one](#id1)
2. [Implementation](#id2)
3. [How to run](#id3)
4. [Understanding the output](#id4)
5. [Intrahost mutations](#id5)
6. [Annotation resources](#id6)
7. [Future work](#id7)
8. [References](#id8)


## Two pipelines in one

In CoVigator we analyse samples from two different sources, ENA and GISAID. While from the first we get the raw reads in
FASTQ format, from the second we obtain already assembled sequences in FASTA format. Each of these formats has to be 
analysed differently. Also, the output data that we can obtain from each of these is different.

### Pipeline for FASTQ files

The FASTQ pipeline includes the following steps:
- **Trimming**. `fastp` is used to trim reads with default values. This step also includes QC filtering.
- **Alignment**. `BWA mem` is used for the alignment of single or paired end samples.
- **BAM preprocessing**. BAM files are prepared and duplicate reads are marked using GATK and Picard tools.
- **Coverage analysis**. `samtools coverage` and `samtools depth` are used to compute the horizontal and vertical 
  coverage respectively.
- **Variant calling**. Four different variant callers are employed: BCFtools, LoFreq, iVar and GATK. 
  Subsequent processing of resulting VCF files is independent for each caller, except for iVar which does not produce a VCF file but a custom TSV file.
- **Variant normalization**. `bcftools norm` and `vt` tools are employed to left align indels, trim variant calls and remove variant duplicates.
- **Variant annotation**. `SnpEff` is employed to annotate the variant consequences of variants, 
  `bcftools annotate` is employed to add additional annotations.

Both single end and paired end FASTQ files are supported.

### Pipeline for FASTA files

The FASTA pipeline includes the following steps:
- **Global Alignment**. A Smith-Waterman global alignment is performed against the reference sequence.
- **Variant calling**. Based on the alignment we call point mutations and small indels. 
  Indels longer than 50 bp and at the beginning or end of the assembly sequence are excluded. Any mutation where
  either reference or assembly contain a N is excluded.
- **Variant normalization**. `bcftools norm` and `vt` tools are employed to left align indels, trim variant calls and remove variant duplicates.
- **Variant annotation**. `SnpEff` is employed to annotate the variant consequences of variants, 
  `bcftools annotate` is employed to add additional annotations.

The FASTA file is expected to contain a single assembly sequence.

## Implementation

The pipeline is implemented as a Nextflow workflow with its DSL2 syntax.
The dependencies are managed through a conda environment to ensure version traceability and reproducibility.
The references for SARS-CoV-2 are embedded in the pipeline.
The pipeline is based on a number of third-party tools, plus a custom implementation based on biopython (Cock, 2009) 
for the alignment and subsequent variant calling over a FASTA file.

All code is open sourced in GitHub [https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline)
and made available under the MIT license. We welcome any contribution. 
If you have troubles using the CoVigator pipeline or you find an issue, we will be thankful if you would report a ticket 
in GitHub.

The alignment, BAM preprocessing and variant normalization pipelines are based on the implementations in additional 
Nextflow pipelines within the TronFlow initiative [https://tronflow-docs.readthedocs.io/](https://tronflow-docs.readthedocs.io/). 


## How to run

### Requirements

- Nextflow >= 20.07.1
- Java >= 8
- Conda >=4.9

### Initializing the conda environments

If you are planning to use it concurrently on multiple samples with conda, first initialize the conda environment, 
otherwise concurrent sample executions will clash trying to create the same conda environment. Run:
```
nextflow main.nf -profile conda --initialize
```

To make sure that all your executions use the same conda environment we recommend using the directive `conda.cacheDir`
[https://www.nextflow.io/docs/latest/config.html#scope-conda](https://www.nextflow.io/docs/latest/config.html#scope-conda).
Otherwise, your conda environment will be stores in the execution folder under `work/conda`.

### Running it

All options available with the `--help` option:
```
$ nextflow run tron-bioinformatics/covigator-ngs-pipeline --help


Usage:
    nextflow run tron-bioinformatics/covigator-ngs-pipeline -profile conda --help

Input:
    * --fastq1: the first input FASTQ file (not compatible with --fasta)
    * --fasta: the FASTA file containing the assembly sequence (not compatible with --fastq1)
    * --name: the sample name, output files will be named after this name
    * --reference: the reference genome FASTA file, *.fai, *.dict and bwa indexes are required.
    * --gff: the GFFv3 gene annotations file (only required with --fastq1)
    * --output: the folder where to publish output
    * --input_fastqs_list: alternative to --name and --fastq1 for batch processing
    * --library: required only when using --input_fastqs
    * --input_fastas_list: alternative to --name and --fasta for batch processing

Optional input:
    * --fastq2: the second input FASTQ file
    * --min_base_quality: minimum base call quality to take a base into account (default: 20)
    * --min_mapping_quality: minimum mapping quality to take a read into account (default: 20)
    * --low_frequency_variant_threshold: VAF threshold to mark a variant as low frequency (default: 0.2)
    * --subclonal_variant_threshold: VAF superior threshold to mark a variant as subclonal (default: 0.8)
    * --memory: the ammount of memory used by each job (default: 3g)
    * --cpus: the number of CPUs used by each job (default: 1)
    * --initialize: initialize the conda environment
    * --skip_lofreq: skips calling variants with LoFreq
    * --skip_gatk: skips calling variants with GATK
    * --skip_bcftools: skips calling variants with BCFTools
    * --skip_ivar: skips calling variants with iVar
    * --match_score: global alignment match score, only applicable for assemblies (default: 2)
    * --mismatch_score: global alignment mismatch score, only applicable for assemblies (default: -1)
    * --open_gap_score: global alignment open gap score, only applicable for assemblies (default: -3)
    * --extend_gap_score: global alignment extend gap score, only applicable for assemblies (default: -0.1)
    * --chromosome: chromosome for variant calls, only applicable for assemblies (default: "MN908947.3")
    * --skip_sarscov2_annotations: skip some of the SARS-CoV-2 specific annotations (default: false)
    * --snpeff_data: path to the SnpEff data folder, it will be useful to use the pipeline on other virus than SARS-CoV-2
    * --snpeff_config: path to the SnpEff config file, it will be useful to use the pipeline on other virus than SARS-CoV-2
    * --snpeff_organism: organism to annotate with SnpEff, it will be useful to use the pipeline on other virus than SARS-CoV-2

Output:
    * Output a normalized, phased and annotated VCF file for each of BCFtools, GATK and LoFreq when FASTQ files are
    provided or a single VCF obtained from a global alignment when a FASTA file is provided
    * Output a TSV file output from iVar
```

It is recommended to specify the pipeline version explicitly for reproducibility purposes:

```
$ nextflow run tron-bioinformatics/covigator-ngs-pipeline -r v0.7.0 ...
```

In order to update your cached version of the pipeline to a specific version use:
```
nextflow pull tron-bioinformatics/covigator-ngs-pipeline -r v0.7.0
```

We recommend using the provided `conda` profile (`-profile conda`), otherwise all dependencies will need to be installed manually and 
made available on the path. In order to combine the `conda` profile with any other custom Nextflow configuration 
(e.g.: a computational cluster queue like Slurm), you will need to use more than one profile. But beware that the 
order matters, the profiles are applied from left to right. A typical situation would be to use the CoVigator conda 
profile plus your default configuration, `-profile conda,standard`.

For paired end reads (ie: two FASTQ files):
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] \
--fastq1 <FASTQ_FILE> \
--fastq2 <FASTQ_FILE> \
--name example_run \
--output <OUTPUT_FOLDER> \
[--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] \
[--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

For single end reads (ie: one FASTQ file):
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] \
--fastq1 <FASTQ_FILE> \
--name example_run \
--output <OUTPUT_FOLDER> \
[--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] \
[--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

For assemblies (ie: FASTA files):
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] \
--fasta <FASTA_FILE> \
--name example_run \
--output <OUTPUT_FOLDER> \
[--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] \
[--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

Sometimes it may be useful to run a batch of samples and not just one by one. The CoVigator pipeline supports this use
case through a tab-separated values file.

For batch processing of reads (ie: FASTQ files) use `--input_fastqs_list` and `--library`: 
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] \
--input_fastqs_list <TSV_FILE> \
--library <paired|single> \
--output <OUTPUT_FOLDER> \
[--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] \
[--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```
where the TSV file contains two or three columns tab-separated columns **without header**. 
Columns: sample name, path to FASTQ 1 and optionally path to FASTQ 2.

| name | FASTQ 1 | FASTQ 2 |
|---------|-----------------|-----------------|
| sample1 | /path/to/sample1.1.fastq | /path/to/sample1.2.fastq |
| sample2 | /path/to/sample2.1.fastq | /path/to/sample2.2.fastq |

Beware that the library only needs to be provided for batch processing. This implies that you cannot process in the 
same batch single and paired-end samples.

For batch processing of assemblies use `--input_fastas_list`.
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] \
--input_fastas_list <TSV_FILE> \
--library <paired|single> \
--output <OUTPUT_FOLDER> \
[--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] \
[--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```
where the TSV file contains two columns tab-separated columns without header. Columns: sample name and path to FASTA.

| name | FASTA |
|---------|-----------------|
| sample1 | /path/to/sample1.fasta |
| sample2 | /path/to/sample2.fasta |


### Using a custom reference

The CoVigator pipeline is applicable to other virus or to other SARS-CoV-2 references. If you want to use CoVigator 
to analyse other virus you will need to specify the reference genome FASTA and GFF files and the chromosome name 
through the parameters `--reference`, `--gff` and `--chromosome`. Notice that it is assumed a single chromosome, 
organisms with multiple chromosomes are not supported.

Additionally, the FASTA needs bwa indexes and .fai index.
These indexes can be generated with the following two commands:
```
bwa index reference.fasta
samtools faidx reference.fasta
```

Furthermore, you will need to provide the adequate SnpEff resources through `--snpeff_data`, `--snpeff_config` and 
`--snpeff_organism`. This is documented in SnpEff 
[https://pcingola.github.io/SnpEff/se_buildingdb/#reference-genome-gtf-gff-refseq-or-genbank](https://pcingola.github.io/SnpEff/se_buildingdb/#reference-genome-gtf-gff-refseq-or-genbank)

Last, but not least, you will need to disable the SARS-CoV-2 specific annotations on conservation and protein domains 
with the flag `--skip_sarscov2_annotations`.

If all this process is obscure to you we will be happy to help you, please, contact us through our GitHub or send us 
an email.

## Understanding the output

Although the VCFs are normalized for both pipelines, the FASTQ pipeline runs four variant callers, while the FASTA
pipeline runs a single variant caller. Also, there are several metrics in the FASTQ pipeline that are not present
in the output of the FASTA pipeline. Here we will describe these outputs.

### FASTQ pipeline output

Find in the table below a description of each of the expected files and a link to a sample file for the FASTQ pipeline.
The VCF files will be described in more detail later.

| Name                        | Description                                                  | Sample file                                                                                          |
|-----------------------------|--------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| $NAME.fastp_stats.json | Output metrics of the fastp trimming process in JSON format | [ERR4145453.fastp_stats.json](_static/covigator_pipeline_sample_output_reads/ERR4145453.fastp_stats.json) |
| $NAME.fastp_stats.html | Output metrics of the fastp trimming process in HTML format | [ERR4145453.fastp_stats.html](_static/covigator_pipeline_sample_output_reads/ERR4145453.fastp_stats.html) |
| $NAME.deduplication_metrics.txt | Deduplication metrics | [ERR4145453.deduplication_metrics.txt](_static/covigator_pipeline_sample_output_reads/ERR4145453.deduplication_metrics.txt) |
| $NAME.coverage.tsv | Coverage metrics (eg: mean depth, % horizontal coverage) | [ERR4145453.coverage.tsv](_static/covigator_pipeline_sample_output_reads/ERR4145453.coverage.tsv) |
| $NAME.depth.tsv | Depth of coverage per position | [ERR4145453.depth.tsv](_static/covigator_pipeline_sample_output_reads/ERR4145453.depth.tsv) |
| $NAME.bcftools.normalized.annotated.vcf.gz | Bgzipped, tabix-indexed and annotated output VCF from BCFtools | [ERR4145453.bcftools.normalized.annotated.vcf.gz](_static/covigator_pipeline_sample_output_reads/ERR4145453.bcftools.normalized.annotated.vcf.gz) |
| $NAME.gatk.normalized.annotated.vcf.gz | Bgzipped, tabix-indexed and annotated output VCF from GATK | [ERR4145453.gatk.normalized.annotated.vcf.gz](_static/covigator_pipeline_sample_output_reads/ERR4145453.gatk.normalized.annotated.vcf.gz) |
| $NAME.lofreq.normalized.annotated.vcf.gz | Bgzipped, tabix-indexed and annotated output VCF from LoFreq | [ERR4145453.lofreq.normalized.annotated.vcf.gz](_static/covigator_pipeline_sample_output_reads/ERR4145453.lofreq.normalized.annotated.vcf.gz) |
| $NAME.ivar.tsv | Output table from iVar | [ERR4145453.ivar.tsv](_static/covigator_pipeline_sample_output_reads/ERR4145453.ivar.tsv) |

**NOTE**: iVar variant caller output variants in a custom format, in particular the indel representation is not 
comparable to the representation of indels in a VCF file. In the future we may parse this format into a VCF file.

### FASTA pipeline output

The FASTA pipeline returns a single VCF file. The VCF files will be described in more detail later.

| Name                        | Description                                                  | Sample file                                                                                          |
|-----------------------------|--------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| $NAME.assembly.normalized.annotated.vcf.gz | Bgzipped, tabix-indexed and annotated output VCF | [ERR4145453.assembly.normalized.annotated.vcf.gz](_static/covigator_pipeline_sample_output_assembly/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz) |

### Understanding the output VCF files

All VCFs from BCFtools, GATK, LoFreq and assemblies are normalized as described in (Tan, 2015). 
Thus their results are comparable. The annotations are in the VCF INFO field.
Although some technical annotations in the VCF are specific to the variant caller, the biological annotations
are shared between all VCFs: ConsHMM conservation scores as reported in (Kwon, 2021), Pfam domains as reported in 
Ensemble annotations and SnpEff effect annotations.

Biological annotations: 

- `INFO/ANN` are the SnpEff consequence annotations (eg: overlapping gene, effect of the mutation). 
This are described in detail here [http://pcingola.github.io/SnpEff/se_inputoutput/](http://pcingola.github.io/SnpEff/se_inputoutput/) 
- `INFO/CONS_HMM_SARS_COV_2` is the ConsHMM conservation score in SARS-CoV-2
- `INFO/CONS_HMM_SARBECOVIRUS` is the ConsHMM conservation score among Sarbecovirus
- `INFO/CONS_HMM_VERTEBRATE_COV` is the ConsHMM conservation score among vertebrate Corona virus
- `INFO/PFAM_NAME` is the Interpro name for the overlapping Pfam domains
- `INFO/PFAM_DESCRIPTION` is the Interpro description for the overlapping Pfam domains

This is an example of biological annotations of a missense mutation in the spike protein on the N-terminal subunit 1 domain.
```
ANN=A|missense_variant|MODERATE|S|gene-GU280_gp02|transcript|TRANSCRIPT_gene-GU280_gp02|protein_coding|1/1|c.118G>A|
p.D40N|118/3822|118/3822|40/1273||;CONS_HMM_SARS_COV_2=0.57215;CONS_HMM_SARBECOVIRUS=0.57215;CONS_HMM_VERTEBRATE_COV=0;
PFAM_NAME=bCoV_S1_N;PFAM_DESCRIPTION=Betacoronavirus-like spike glycoprotein S1, N-terminal
```

The different technical annotations from the different variant callers are out of the scope of this documentation.
For illustration purpose see an example from each of the callers on the FASTQ pipeline.

- BCFtools: `DP=137;VDB=0.648575;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,109,23;MQ=60;`
- GATK: `AC=1;AF=1;AN=1;BaseQRankSum=2.481;DP=261;FS=2.061;MLEAC=1;MLEAF=1;MQ=60;MQRankSum=0;QD=28.2;ReadPosRankSum=2.858;SOR=0.23;`
- LoFreq: `DP=245;AF=0.020408;SB=0;DP4=113,127,2,3`


## Intrahost mutations

Some mutations may be observed in a subset of the virus sample, this may arise through intrahost virus evolution or
co-infection. These mutations observed at a frequency lower than 100 % can only be detected when analysing the raw 
reads (ie: the FASTQs) as in the assembly (ie: the FASTA file) a single virus consensus sequence is represented.

These intrahost lower variant allele frequency (VAF) mutations are challenging to distinguish from sequencing and
analytical errors. The variant caller best suited for this purpose is LoFreq. 
Also, LoFreq provides among its technical annotations the VAF and depth of coverage of every mutation.
- `INFO/DP` is the number of reads overlapping this position
- `INFO/AF` is the VAF as reported by LoFreq (only avalable in LoFreq calls)

The LoFreq variants are annotated on the `FILTER` column using the reported VAF into three categories: 
- `LOW_FREQUENCY`
- `SUBCLONAL`
- `PASS` variants corresponding to clonal variants 

By default, variants with a VAF < 20 % are considered `LOW_FREQUENCY` and variants with a VAF >= 20 % and < 80 % are considered 
`SUBCLONAL`. This thresholds can be changed with the parameters `--low_frequency_variant_threshold` and
`--subclonal_variant_threshold`.

Thus all variants with a VAF lower than 0.8 could be considered as intrahost variants. Beware that the lower the VAF
the higher the enrichment for false positive calls.

## Annotations resources

SARS-CoV-2 ASM985889v3 references were downloaded from Ensembl on 6th of October 2020:
- ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
- ftp://ftp.ensemblgenomes.org/pub/viruses/gff3/sars_cov_2/Sars_cov_2.ASM985889v3.101.gff3.gz

ConsHMM mutation depletion scores downloaded on 1st of July 2021:
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionConsHMM.bed
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionSarbecovirusConsHMM.bed
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionVertebrateCoVConsHMM.bed

Gene annotations including Pfam domains downloaded from Ensembl on 25th of February 2021 from:
- ftp://ftp.ensemblgenomes.org/pub/viruses/json/sars_cov_2/sars_cov_2.json


## Future work

- Normalization of iVar results into a standard annotated VCF file
- Validation of intrahost variant calls
- Pipeline for Oxford Nanopore technology
- Homogeneisation of the technical annotations between variant callers
- Annotation of variants with known lineages + pangolin integration
- Variant calls from assemblies contain an abnormally high number of deletions of size greater than 3 bp. This
is a technical artifact that would need to be avoided


## References

- Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. Unified Representation of Genetic Variants. Bioinformatics (2015) 31(13): 2202-2204](http://bioinformatics.oxfordjournals.org/content/31/13/2202) and uses bcftools [Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (Oxford, England), 27(21), 2987–2993. 10.1093/bioinformatics/btr509
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.
- Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M. (2013). From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. Curr Protoc Bioinformatics, 43:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.
- Martin, M., Patterson, M., Garg, S., O Fischer, S., Pisanti, N., Klau, G., Schöenhuth, A., & Marschall, T. (2016). WhatsHap: fast and accurate read-based phasing. BioRxiv, 085050. https://doi.org/10.1101/085050
- Danecek, P., & McCarthy, S. A. (2017). BCFtools/csq: haplotype-aware variant consequences. Bioinformatics, 33(13), 2037–2039. https://doi.org/10.1093/bioinformatics/btx100
- Wilm, A., Aw, P. P. K., Bertrand, D., Yeo, G. H. T., Ong, S. H., Wong, C. H., Khor, C. C., Petric, R., Hibberd, M. L., & Nagarajan, N. (2012). LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Research, 40(22), 11189–11201. https://doi.org/10.1093/nar/gks918
- Grubaugh, N. D., Gangavarapu, K., Quick, J., Matteson, N. L., De Jesus, J. G., Main, B. J., Tan, A. L., Paul, L. M., Brackney, D. E., Grewal, S., Gurfield, N., Van Rompay, K. K. A., Isern, S., Michael, S. F., Coffey, L. L., Loman, N. J., & Andersen, K. G. (2019). An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biology, 20(1), 8. https://doi.org/10.1186/s13059-018-1618-7
- Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
- Kwon, S. Bin, & Ernst, J. (2021). Single-nucleotide conservation state annotation of the SARS-CoV-2 genome. Communications Biology, 4(1), 1–11. https://doi.org/10.1038/s42003-021-02231-w
- Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.
