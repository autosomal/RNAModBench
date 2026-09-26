# How the tool inventory was assembled

`tools/inventory/tables/` was built to answer a specific question: for each tool,
what exactly was run, with what version, model, thresholds and coordinate
handling. Three evidence layers were combined, and the tables record which layer
each value came from.

## Layers

1. **Automatic collection** (`tools/inventory/scripts/`, re-runnable):
   * `collect_commands.py` scans every shell driver that was used for a run and
     extracts the command lines with their file and line number → `TI2`.
   * `collect_params.py` parses the YAML/environment specifications, run logs and
     notes for parameters and model references.
   * `collect_versions.py` probes versions three ways — `conda list` per
     environment, each binary's `--version`, and the source checkout — → `TI3`, `TI4`.
2. **Manual curation** — `tools/inventory/curated/tool_inventory_curated.csv` is the
   only file that was edited by hand; the tables are generated from it
   (`scripts/build_tables.py`).
3. **Per-tool reading** — where a value was not present in a log or a YAML it was
   read from the tool's own argparse/click defaults, README or model file, and the
   reading was recorded with its evidence pointer
   (`scripts/parse_research.py` → `apply_research.py`, promoted to `CONFIRMED` by
   `promote_status.py` after human review).

## What the recorded command lines are, and are not

The `method` column of `TI2_command_lines.csv` is `static_script` for all 3,608
rows: the command lines were **extracted from the scripts that ran them**, not
captured from a process log. That has three consequences worth stating:

* a line shows the invocation as written, so shell variables appear as the script
  spelled them rather than as their expanded values;
* where a driver looped over samples, one line represents the loop body and the
  `dataset` column says which iteration it describes;
* lines that only echo, log or test are included, because dropping them would have
  required interpreting the script.

Paths inside the commands point at the working tree of the machine that ran them.
They have been rewritten to `$RNAMODBENCH_ROOT/...` (this repository) or
`$RNAMODBENCH_LOCAL/...` (inputs that are not redistributed); the flags, arguments
and ordering are unchanged. Legacy pre-reorganisation prefixes such as
`result/` and `code_user/` are kept as they were written, so an old line still
matches the script it came from.

## Versions and models

`TI3_software_versions.csv` carries the probing method per row
(`conda_list`, `bin_version`, `static_yml`, `source_git`, …), so a value can be
weighed by how it was obtained. `TI4_model_checkpoints.csv` names the model or
checkpoint file each tool configuration used. `conda_env` is blank on the rows
where the running environment could not be established from the script.

## Coverage and its limits

`TI1_per_tool_implementation.csv`/`.md` give the ten requested fields per tool
configuration, each with an evidence pointer, and `TI8_coverage_report.md` is the
completeness matrix. Two limits are inherent rather than fixable here:

* **Configuration files exist for only two tools.** xPore and CHEUI-diff were
  driven from YAML files, and those files are deposited under `tools/configs/`
  with their paths rewritten the same way. Everything else was driven by
  command-line flags, which is what `TI2` records.
* **Not every conda environment has a committed specification.** 11 tool
  environment specifications are in `envs/as_run/` (and the pipeline's own `envs/*.yaml`); the run scripts
  reference more environments than that, and for those the version evidence is
  the `conda_list` rows in `TI3` rather than a lock file.

`tools/inventory/scripts/run_all.sh` re-runs the collection end to end; it needs
the tool installations and the result trees under `$RNAMODBENCH_LOCAL`, so it is a
maintenance script for the authors rather than part of a figure rebuild.
