from argparse import ArgumentParser

from covigator.precomputations.load_cooccurrences import CooccurrenceMatrixLoader
from covigator.precomputations.loader import PrecomputationsLoader
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
import covigator
import covigator.configuration
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.accessor.gisaid_accessor import GisaidAccessor
from covigator.configuration import Configuration
from covigator.database.database import Database
from covigator.pipeline.ena_pipeline import Pipeline
from covigator.pipeline.gisaid_pipeline import GisaidPipeline
from covigator.processor.ena_downloader import EnaDownloader
from covigator.processor.ena_processor import EnaProcessor
from covigator.processor.gisaid_processor import GisaidProcessor
from logzero import logger


def ena_accessor():
    parser = ArgumentParser(
        description="Covigator {} ENA accessor".format(covigator.VERSION))
    parser.add_argument(
        "--tax-id",
        dest="tax_id",
        help="the taxonomy id of the organism to analyse, eg: 2697049 for Sars-COV-2",
        default="2697049"
    )
    parser.add_argument(
        "--host-tax-id",
        dest="host_tax_id",
        help="the taxonomy id of the host organism, eg: 9606 for Homo sapiens",
        default="9606"
    )

    args = parser.parse_args()
    tax_id = args.tax_id
    host_tax_id = args.host_tax_id
    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_accesor)
    EnaAccessor(tax_id=tax_id, host_tax_id=host_tax_id, database=Database(config=config, initialize=True)).access()


def gisaid_accessor():
    parser = ArgumentParser(
        description="Covigator {} GISAID accessor".format(covigator.VERSION))
    parser.add_argument(
        "--input-fasta",
        dest="input_fasta",
        help="the FASTA input file from GISAID",
        required=True
    )
    parser.add_argument(
        "--input-metadata",
        dest="input_metadata",
        help="GISAID metadata file",
        required=True
    )

    args = parser.parse_args()

    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_accesor)
    GisaidAccessor(input_fasta=args.input_fasta, input_metadata=args.input_metadata, database=Database(initialize=True, config=config)).access()


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
        "--num-jobs",
        dest="num_jobs",
        help="The number of dask jobs to spin, this corresponds to the number of whole nodes requested to the cluster. "
             "The configuration of each job is in the correspoding jobqueue.yaml file",
        default=1
    )
    parser.add_argument(
        "--local",
        dest="local",
        help="if set it runs locally, no slurm cluster",
        action='store_true',
        default=False
    )
    parser.add_argument(
        "--num-local-cpus",
        dest="num_local_cpus",
        help="number of CPUs to be used by the processor when running locally",
        default=1
    )
    parser.add_argument(
        "--download",
        dest="download",
        help="if set it tries to download ENA FASTQs if they are not already in place. Not applicale for GISAID",
        action='store_true',
        default=False
    )

    args = parser.parse_args()
    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_processor)
    if args.local:
        _start_dask_processor(args, config, args.download, num_local_cpus=int(args.num_local_cpus))
    else:
        with SLURMCluster(
                walltime='72:00:00',  # hard codes maximum time to 72 hours
                scheduler_options={"dashboard_address": ':{}'.format(config.dask_port)}) as cluster:
            cluster.scale(jobs=int(args.num_jobs))
            _start_dask_processor(args, config, args.download, cluster=cluster)


def ena_downloader():
    parser = ArgumentParser(
        description="Covigator {} ENA downloader".format(covigator.VERSION))

    config = Configuration()
    covigator.configuration.initialise_logs(config.logfile_accesor)
    EnaDownloader(database=Database(config=config, initialize=True), config=config).process()


def _start_dask_processor(args, config, download, cluster=None, num_local_cpus=1):
    with Client(cluster) if cluster is not None else Client(n_workers=num_local_cpus, threads_per_worker=1) as client:
        if args.data_source == "ENA":
            EnaProcessor(database=Database(initialize=True, config=config),
                         dask_client=client, config=config, download=download) \
                .process()
        elif args.data_source == "GISAID":
            GisaidProcessor(database=Database(initialize=True, config=config),
                            dask_client=client, config=config) \
                .process()
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
    vcf_file, qc_file = Pipeline(config=Configuration()).run(
        run_accession="test", fastq1=args.fastq1, fastq2=args.fastq2)
    logger.info("Output VCF file: {}".format(vcf_file))
    logger.info("Output QC file: {}".format(qc_file))


def gisaid_pipeline():
    parser = ArgumentParser(description="Run Pipeline for testing")
    parser.add_argument("--run_accession", dest="run_accession", help="Specify run accession to process", required=True)

    args = parser.parse_args()
    # FIXME: this interface is broken, it requires the SampleGisaid object not the run_accession
    vcf_file = GisaidPipeline(config=Configuration()).run(run_accession=args.run_accession)
    logger.info("Output VCF file: {}".format(vcf_file))


def precompute_queries():
    parser = ArgumentParser(description="Precompute some aggregation queries")
    parser.parse_args()

    database = Database(initialize=True, config=Configuration())
    loader = PrecomputationsLoader(session=database.get_database_session())
    logger.info("Starting precomputation...")
    loader.load()
    logger.info("Done precomputing")


def cooccurrence():
    parser = ArgumentParser(description="Precompute cooccurrence of mutations")
    parser.add_argument(
        "--source",
        dest="data_source",
        help="Specify data source. This can be either ENA or GISAID",
        required=True
    )
    parser.add_argument(
        "--maximum-mutation-length",
        dest="maximum_length",
        help="Only mutations with this maximum size will be included in the cooccurence matrix",
        default=10
    )
    args = parser.parse_args()

    database = Database(initialize=True, config=Configuration())
    loader = CooccurrenceMatrixLoader(session=database.get_database_session(), source=args.data_source)
    logger.info("Starting precomputation...")
    loader.load(maximum_length=int(args.maximum_length))
    logger.info("Done precomputing")
