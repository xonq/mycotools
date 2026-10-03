#! /usr/bin/env python3
"""Modify a linked primary MycotoolsDB's persistent configuration.

The options that shape what future `mtdb update` runs acquire -- inclusion of
use-restricted data, MycoCosm (JGI) participation, and lineage constraints --
are established when the database is initialized (`mtdb update -i`). Because
those settings define the database, `mtdb update` only accepts them at
initialization; this utility is the supported way to change them afterward.

It edits `config/mtdb.json` in place. With no options it prints the current
configuration. Every change takes effect on the next `mtdb update`."""

import os
import sys
import logging
import argparse
from collections import defaultdict
from pathlib import Path
from mycotools.lib.kontools import (
    format_path,
    read_json,
    write_json,
    split_input,
    setup_logging,
)

logger = logging.getLogger(__name__)

# ranks accepted for lineage constraints; mirrors `mtdb update`
PERMITTED_RANKS = {"phylum", "subphylum", "class", "order", "family", "genus"}

# config keys written by update.gen_config, with a label for display. Ordered
# most- to least- frequently adjusted.
CONFIG_LABELS = (
    ("branch", "Kingdom/branch"),
    ("jgi", "MycoCosm (JGI)"),
    ("nonpublished", "Use-restricted data"),
    ("lineage_constraints", "Lineage constraints"),
    ("rogue", "Standalone (rogue)"),
    ("repository", "Reference repository"),
    ("forbidden", "Forbidden-ome ledger"),
)

# values `nonpublished` may take when enabled (config stores "yes"; the historic
# forms are accepted defensively)
_TRUE = {"yes", "y", "true"}


def load_config():
    """Return (config_path, config) for the linked primary MTDB, or exit.

    The configuration lives beside the linked database, so a missing MYCODB is
    "nothing is linked" and a missing file is a corrupt install -- the two
    exits mirror the codes `mtdb update` raises for the same conditions."""
    if "MYCODB" not in os.environ:
        logger.error("MTDB not linked. Link via `mtdb -i <DB_PATH>`")
        sys.exit(50)
    path = format_path("$MYCODB/../config/mtdb.json")
    if not Path(path).is_file():
        logger.error("corrupted MycotoolsDB - no configuration found")
        sys.exit(21)
    return path, read_json(path)


def _is_true(value):
    """Whether a stored `nonpublished`/`jgi` value counts as enabled."""
    if isinstance(value, str):
        return value.lower() in _TRUE
    return bool(value)


def fmt_value(key, value):
    """Render a config value for display."""
    if key == "lineage_constraints":
        if not value:
            return "none"
        return "; ".join(
            f"{rank}: {', '.join(sorted(lineages))}"
            for rank, lineages in sorted(value.items())
        )
    if key in {"nonpublished", "jgi", "rogue"}:
        return "yes" if _is_true(value) else "no"
    if value in (None, ""):
        return "none"
    return str(value)


def show_config(config):
    """Print the current configuration."""
    print("MycotoolsDB configuration:", flush=True)
    for key, label in CONFIG_LABELS:
        if key in config:
            print(f"  {label}: {fmt_value(key, config[key])}", flush=True)


def parse_lineage_constraints(lineage, rank):
    """Turn positional --lineage/--rank strings into a rank->[lineages] dict.

    Mirrors the parsing in `mtdb update`'s control_flow so a constraint set
    reads identically whether it is established at init or reconfigured here.
    Exits on a length mismatch or an unrecognized rank."""
    lineage_constraints = split_input(lineage)
    rank_constraints = split_input(rank)
    if len(lineage_constraints) != len(rank_constraints):
        logger.error("--lineage must be same length as --rank")
        sys.exit(18)
    rank2lineages = defaultdict(set)
    for i, lin in enumerate(lineage_constraints):
        rank_c = rank_constraints[i].lower()
        if rank_c not in PERMITTED_RANKS:
            logger.error(f"accepted ranks: {sorted(PERMITTED_RANKS)}")
            sys.exit(22)
        rank2lineages[rank_c].add(lin.lower())
    return {k: sorted(v) for k, v in sorted(rank2lineages.items())}


def _set_nonpublished(config, args, changes):
    """Toggle use-restricted data inclusion in `config`."""
    if args.nonpublished and args.published:
        logger.error("--nonpublished and --published are mutually exclusive")
        sys.exit(1)
    if args.nonpublished:
        if _is_true(config.get("nonpublished")):
            logger.info("use-restricted data already enabled")
            return
        # the T&C acknowledgement lives in update; import it lazily so display
        # and the other toggles do not pull in the download stack
        from mycotools.mtdb.update import validate_t_and_c

        config["nonpublished"] = validate_t_and_c(config, discrepancy=True)
        changes.append("use-restricted data enabled")
    elif args.published:
        if not _is_true(config.get("nonpublished")):
            logger.info("use-restricted data already disabled")
            return
        config["nonpublished"] = False
        changes.append("use-restricted data disabled")


def _set_jgi(config, args, changes):
    """Toggle MycoCosm (JGI) participation in `config`."""
    if args.ncbi_only and args.jgi:
        logger.error("--ncbi_only and --jgi are mutually exclusive")
        sys.exit(1)
    if args.ncbi_only:
        if not _is_true(config.get("jgi")):
            logger.info("MycoCosm already disabled")
            return
        config["jgi"] = False
        changes.append("MycoCosm disabled (NCBI only)")
    elif args.jgi:
        if config.get("branch", "fungi") != "fungi":
            logger.warning(
                f'branch "{config.get("branch")}" does not use MycoCosm; '
                "enabling it has no effect"
            )
        if _is_true(config.get("jgi")):
            logger.info("MycoCosm already enabled")
            return
        config["jgi"] = True
        changes.append("MycoCosm enabled")


def _set_lineage(config, args, changes):
    """Set or clear the lineage constraints in `config`."""
    if args.clear_lineage and (args.lineage or args.rank):
        logger.error("--clear_lineage cannot be combined with --lineage/--rank")
        sys.exit(1)
    if args.clear_lineage:
        if config.get("lineage_constraints"):
            config["lineage_constraints"] = {}
            changes.append("lineage constraints cleared")
        else:
            logger.info("no lineage constraints to clear")
    elif args.lineage or args.rank:
        if not (args.lineage and args.rank):
            logger.error("--lineage requires --rank")
            sys.exit(16)
        new = parse_lineage_constraints(args.lineage, args.rank)
        if new != config.get("lineage_constraints"):
            config["lineage_constraints"] = new
            changes.append("lineage constraints updated")
        else:
            logger.info("lineage constraints unchanged")


def apply_changes(config, args):
    """Mutate `config` in place per `args`; return a list of change summaries."""
    changes = []
    _set_nonpublished(config, args, changes)
    _set_jgi(config, args, changes)
    _set_lineage(config, args, changes)
    return changes


def cli():
    parser = argparse.ArgumentParser(
        description="Modify the linked primary MycotoolsDB configuration "
        "(config/mtdb.json). With no options, print the current configuration. "
        "Changes take effect on the next `mtdb update`."
    )
    restr = parser.add_argument_group("Use-restricted data")
    restr.add_argument(
        "--nonpublished",
        action="store_true",
        help="[FUNGI]: Include MycoCosm use-restricted data",
    )
    restr.add_argument(
        "--published",
        action="store_true",
        help="Exclude use-restricted data",
    )

    jgi = parser.add_argument_group("MycoCosm")
    jgi.add_argument(
        "--ncbi_only", action="store_true", help="[FUNGI]: Forego MycoCosm (NCBI only)"
    )
    jgi.add_argument("--jgi", action="store_true", help="[FUNGI]: Include MycoCosm")

    lin = parser.add_argument_group("Lineage constraints")
    lin.add_argument("-l", "--lineage", help="Lineage(s) to constrain updates to")
    lin.add_argument(
        "-rk", "--rank", help="Rank(s) that positionally correspond to -l"
    )
    lin.add_argument(
        "--clear_lineage", action="store_true", help="Remove all lineage constraints"
    )
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    path, config = load_config()

    requested = (
        args.nonpublished
        or args.published
        or args.ncbi_only
        or args.jgi
        or args.lineage
        or args.rank
        or args.clear_lineage
    )
    if not requested:
        show_config(config)
        sys.exit(0)

    changes = apply_changes(config, args)
    if changes:
        write_json(config, path)
        for change in changes:
            logger.info(change)
        logger.info("Run `mtdb update` to apply")
    else:
        logger.info("No configuration changes")
    show_config(config)
    sys.exit(0)


def main():
    cli()


if __name__ == "__main__":
    cli()
