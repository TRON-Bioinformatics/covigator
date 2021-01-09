from argparse import ArgumentParser

from dask.distributed import Client

import covigator
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.dashboard.dashboard import app
from covigator.model import Database
from covigator.processor.processor import Processor


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
    EnaAccessor(tax_id=tax_id, host_tax_id=host_tax_id, database=Database()).access()


def processor():
    parser = ArgumentParser(
        description="Covigator {} processor".format(covigator.VERSION))
    parser.add_argument(
        "--num-cpus",
        dest="num_cpus",
        help="number of CPUs to be used by the processor, this reflects the number of jobs that will be processed in "
             "parallel",
        required=True
    )

    args = parser.parse_args()

    client = Client(n_workers=args.num_cpus, threads_per_worker=1)
    Processor(database=Database(), dask_client=client).process()


def dashboard():
    app.run_server(debug=True)
    #Dashboard(database=Database()).run(debug=True)
