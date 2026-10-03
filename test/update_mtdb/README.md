# `mtdb update` argument-scenario tests

Fast, offline tests for `mycotools/mtdb/update.py` (the `mtdb update` / `mtdb u`
subcommand). They exercise the layers a user drives through the command line,
and deliberately stop at the JGI/NCBI login boundary — **no credentials,
downloads, or genome assimilation.**

## Run

```bash
micromamba run -n mycotools pytest test/update_mtdb/ -v
```

(~2.5 s, 43 tests.)

## What is covered

| Area | How | Examples |
| --- | --- | --- |
| `control_flow` argument validation | called directly with a dummy `ncbi_email` (which short-circuits `loginCheck`); every guard fires before login | bad flag combos → documented exit codes (13–22, 431) |
| Valid/invalid `--kingdom` | direct calls | abbreviations + full names accepted; junk → 431 |
| argparse CLI surface | `python -m mycotools.mtdb.update` subprocess (parser exits before the `datasets` check) | `--help`, unknown flag, non-int `--cpu`/`--resume` |
| Login-free helpers | driven with `test/ust.mtdb` + an in-code genus lineage filter | `gen_config`, `internal_redundancy_check` dereplication, genus-only `extract_constraint_lineages`, ledger/sidecar round-trips |

## What is intentionally **not** covered

Anything past the login boundary: `ref_update`, `rogue_update`,
`taxonomy_update`, MycoCosm/NCBI downloads, Entrez taxonomy queries, and the
`datasets` utility. These need authenticated network access and would download
whole genomes.

## Notes for maintainers

- The validation exit codes are `mtdb update`'s own contract; a change there is
  a behavior change and should update `VALIDATION_CASES` deliberately.
- `--predb` alone exits `15` (not `18`) because the first guard does not count
  `--predb` as a mode — this is asserted (`predb-alone`) to lock in the current
  behavior.
- `MYCOFNA`/`MYCOFAA`/`MYCOGFF3` are set to dummy prefixes by an autouse fixture;
  the pandas MTDB helpers only string-concatenate them, so no real paths exist.
