import os
from argparse import ArgumentParser
from dask.distributed import Client
import covigator
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.database.database import Database
from covigator.processor.pipeline import Pipeline
from covigator.processor.gisaid_pipeline import GisaidPipeline
from covigator.processor.ena_processor import EnaProcessor
from covigator.processor.gisaid_processor import GisaidProcessor
import logzero
from logzero import logger


def ena_accessor():
    parser = ArgumentParser(
        description="Covigator {} ENA accessor".format(covigator.VERSION))
    parser.add_argument(
        "--tax-id",
        dest="tax_id",
        help="the taxonomy id of the organism to analyse, eg: 2697049 for Sars-COV-2",
        required=True
    )
    parser.add_argument(
        "--host-tax-id",
        dest="host_tax_id",
        help="the taxonomy id of the host organism, eg: 9606 for Homo sapiens",
        required=True
    )

    args = parser.parse_args()
    tax_id = args.tax_id
    host_tax_id = args.host_tax_id
    log_file = os.getenv(covigator.ENV_COVIGATOR_ACCESSOR_LOG_FILE)
    if log_file is not None:
        logzero.logfile(log_file, maxBytes=1e6, backupCount=3)
    EnaAccessor(tax_id=tax_id, host_tax_id=host_tax_id, database=Database()).access()


def processor():
    parser = ArgumentParser(
        description="Covigator {} processor".format(covigator.VERSION))
    parser.add_argument(
        "--source",
        dest="data_source",
        help="Specify data source. This can be either ENA or GISAID",
        required=True
    )
    parser.add_argument(
        "--num-cpus",
        dest="num_cpus",
        help="number of CPUs to be used by the processor, this reflects the number of jobs that will be processed in "
             "parallel",
        required=True
    )

    args = parser.parse_args()
    log_file = os.getenv(covigator.ENV_COVIGATOR_PROCESSOR_LOG_FILE)
    if log_file is not None:
        logzero.logfile(log_file, maxBytes=1e6, backupCount=3)
    client = Client(n_workers=int(args.num_cpus), threads_per_worker=1)
    if args.data_source == "ENA":
        EnaProcessor(database=Database(), dask_client=client).process()
    elif args.data_source == "GISAID":
        GisaidProcessor(database=Database(), dask_client=client).process()
    else:
        logger.error("Unknown data source. Please choose either ENA or GISAID")


def pipeline():
    parser = ArgumentParser(description="Run Pipeline for testing")
    parser.add_argument("--fastq1", dest="fastq1",
                        help="First FASTQ to process. In case of single end sequencing this will be the only FASTQ,"
                             "in the case of paired end sequencing this will be the first from the pair", required=True)
    parser.add_argument("--fastq2", dest="fastq2",
                        help="Second FASTQ to process for paired end sequencing, otherwise leave empty")

    args = parser.parse_args()
    vcf_file = Pipeline().run(fastq1=args.fastq1, fastq2=args.fastq2)
    logger.info("Output VCF file: {}".format(vcf_file))


def gisaid_pipeline():
    parser = ArgumentParser(description="Run Pipeline for testing")
    parser.add_argument("--run_accession", dest="run_accession", help="Specify run accession to process", required=True)

    args = parser.parse_args()
    vcf_file = GisaidPipeline().run(run_accession=args.run_accession)
    logger.info("Output VCF file: {}".format(vcf_file))
