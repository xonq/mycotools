#! /usr/bin/env python3
"""Shared factory for nested subcommand dispatchers.

Mycotools groups related tools under a parent command (e.g. `mycotools
download`, `mtdb accession`). A group dispatcher forwards
`<prog> <SUBCOMMAND> ...` to a submodule that parses its own `sys.argv`.

Rather than reimplement the same forward-argv-and-delegate logic in every
group's `__init__.py`, each group builds a `Dispatcher` with its subcommand map
and help text, then exposes the dispatcher's `main`/`cli`. The one exception is
the `mtdb` root command, which carries extra base operations (link/unlink/list,
ome lookup) and stays bespoke."""
import sys
import logging
import argparse
import importlib

logger = logging.getLogger(__name__)


class Dispatcher:
    """A nested subcommand dispatcher for a single command group.

    Parameters
    ----------
    prog : str
        Program name shown in help and prepended to the forwarded argv, e.g.
        ``"mycotools download"``.
    package : str
        Importable package the submodules live in, e.g.
        ``"mycotools.download"``. A resolved subcommand ``jgi`` is run by
        importing ``mycotools.download.jgi`` and calling its ``cli()``.
    subcommands : dict
        Map of subcommand name/alias -> submodule name within ``package``.
    description : str
        Full ``argparse`` description (RawDescriptionHelpFormatter), typically
        a header plus a formatted subcommand listing.
    metavar : str
        Positional metavar shown in usage (default ``"SUBCOMMAND"``).
    arg_help : str
        Help string for the positional (default ``"subcommand to run (see below)"``).
    """

    def __init__(
        self,
        prog,
        package,
        subcommands,
        description,
        metavar="SUBCOMMAND",
        arg_help="subcommand to run (see below)",
    ):
        self.prog = prog
        self.package = package
        self.subcommands = subcommands
        self.description = description
        self.metavar = metavar
        self.arg_help = arg_help

    def build_parser(self):
        """Build the group parser.

        The subcommand positional is followed by a REMAINDER so the submodule's
        arguments (including forwarded flags such as `-h`) pass through verbatim;
        native argparse subparsers would intercept them."""
        parser = argparse.ArgumentParser(
            prog=self.prog,
            description=self.description,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument(
            "subcommand", nargs="?", metavar=self.metavar, help=self.arg_help
        )
        parser.add_argument("rest", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
        return parser

    def delegate(self, module_name, args):
        """Run a submodule's CLI in-process and return its exit code.

        The submodule parses `sys.argv` itself, so its argv is swapped in for the
        call and the `SystemExit` it raises (argparse errors, explicit exits) is
        translated back to an exit code. Imports are deferred so a bare group
        invocation stays light."""
        module = importlib.import_module(f"{self.package}.{module_name}")
        saved_argv = sys.argv
        sys.argv = [f"{self.prog} {module_name}"] + list(args)
        try:
            module.cli()
            return 0
        except SystemExit as exc:
            if exc.code is None:
                return 0
            return exc.code if isinstance(exc.code, int) else 1
        finally:
            sys.argv = saved_argv

    def main(self, argv=None):
        if argv is None:
            argv = sys.argv
        parser = self.build_parser()
        args = parser.parse_args(argv[1:])
        if args.subcommand in self.subcommands:
            sys.exit(self.delegate(self.subcommands[args.subcommand], args.rest))
        if args.subcommand is not None:
            logger.error(f"invalid subcommand: {args.subcommand}")
        parser.print_help(sys.stderr)
        sys.exit(1)

    def cli(self):
        self.main(sys.argv)
