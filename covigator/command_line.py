from argparse import ArgumentParser
import covigator
from covigator.accessor.ena_accessor import EnaAccessor
from covigator.model import Database


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
    accessor = EnaAccessor(tax_id=tax_id, host_tax_id=host_tax_id, database=Database())
    accessor.access()


def pipeline():
    raise NotImplementedError()
