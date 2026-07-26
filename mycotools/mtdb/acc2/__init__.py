#! /usr/bin/env python3
"""Accession-retrieval dispatcher for the `mtdb accession` subcommand.

Routes `mtdb accession <FORMAT> ...` (equivalently `mtdb a <FORMAT> ...`) to a
per-format retrieval module in this package. Each format module parses its own
`sys.argv`, so its argv is swapped in and run in-process; deferred imports keep a
bare `mtdb accession` light."""
import sys
import logging
import argparse
import importlib

logger = logging.getLogger(__name__)

# format name/alias -> submodule within this package (mycotools.mtdb.acc2.<module>)
SUBCOMMANDS = {
    "fa": "fa",
    "fasta": "fa",
    "gff": "gff",
    "gff3": "gff",
    "gbk": "gbk",
    "genbank": "gbk",
    "locus": "locus",
    "l": "locus",
}

DESCRIPTION = """Retrieve MTDB data for accession(s) in a given file format

Formats (all following arguments are forwarded to the format's retriever):
  fa     (fasta)     retrieve protein/nucleotide FASTA for accession(s)
  gff    (gff3)      retrieve GFF3 entries for accession(s)
  gbk    (genbank)   generate a GenBank file for accession(s)/ome(s)
  locus  (l)         retrieve the locus surrounding accession(s)

Examples:
  mtdb accession fa -a <OME>_<ACC>
  mtdb accession gff -h"""


def build_parser():
    """Build the `mtdb accession` argument parser.

    The format positional is followed by a REMAINDER so the format module's
    arguments pass through verbatim (native argparse subparsers would intercept
    forwarded flags such as `-h`)."""
    parser = argparse.ArgumentParser(
        prog="mtdb accession",
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "format",
        nargs="?",
        metavar="FORMAT",
        help="retrieval format (see below)",
    )
    parser.add_argument("rest", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
    return parser


def delegate(module_name, args):
    """Run a format retriever's CLI in-process and return its exit code.

    The format modules live in this package (`mycotools.mtdb.acc2.<module>`) and
    parse `sys.argv` themselves, so their argv is swapped in for the call and the
    `SystemExit` they raise (argparse errors, explicit exits) is translated back
    to an exit code."""
    module = importlib.import_module(f"mycotools.mtdb.acc2.{module_name}")
    saved_argv = sys.argv
    sys.argv = [f"mtdb accession {module_name}"] + list(args)
    try:
        module.cli()
        return 0
    except SystemExit as exc:
        if exc.code is None:
            return 0
        return exc.code if isinstance(exc.code, int) else 1
    finally:
        sys.argv = saved_argv


def main(argv=sys.argv):
    parser = build_parser()
    args = parser.parse_args(argv[1:])
    if args.format in SUBCOMMANDS:
        sys.exit(delegate(SUBCOMMANDS[args.format], args.rest))
    if args.format is not None:
        logger.error(f"invalid format: {args.format}")
    parser.print_help(sys.stderr)
    sys.exit(1)


def cli():
    main(sys.argv)


if __name__ == "__main__":
    cli()
