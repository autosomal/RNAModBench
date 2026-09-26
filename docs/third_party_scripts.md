# Tool scripts the Snakefile expects but does not ship

The pipeline rules invoke helper scripts that belong to the detection tools
themselves — they are the tools' own code, under the tools' own licences, so they
are **not** copied into this repository. The workflow files reference them under
`scripts/`; obtain the tool, then place (or symlink) the file under the name the
rule expects.

| name the Snakefile uses | ships with | file in the tool's distribution |
|---|---|---|
| `CHEUI_preprocess_m6A.py` | CHEUI ([github.com/comprna/CHEUI](https://github.com/comprna/CHEUI)) | `scripts/CHEUI_preprocess_m6A.py` |
| `CHEUI_predict_model1.py` | CHEUI | `scripts/CHEUI_predict_model1.py` |
| `CHEUI_predict_model2.py` | CHEUI | `scripts/CHEUI_predict_model2.py` |
| `Epinano_Variants.py` | EpiNano 1.2.0 ([github.com/enovoa/EpiNano](https://github.com/enovoa/EpiNano)) | `Epinano_Variants.py` |
| `Epinano_Predict.py` | EpiNano 1.2.0 | `Epinano_Predict.py` |
| `Slide_Variants.py` | EpiNano 1.2.0 | `misc/Slide_Variants.py` |
| `DENA_extract.py` | DENA | `step4_predict/LSTM_extract.py` |
| `DENA_LSTM_predict.py` | DENA | `step4_predict/LSTM_predict.py` |
| `MINES_cDNA.py` | MINES | `cDNA_MINES1.py` |
| `extract_raw_and_feature_fast_AUCG.py` | NanoNm | `extract_raw_and_feature_fast_AUCG.py` |
| `predict_sites_Nm.final.py` | NanoNm | `predict_sites_Nm.final.py` |

Two of these names differ from the upstream file name (DENA and MINES); the table
gives the mapping that was actually used in the reported runs.

The scripts that *are* in `scripts/` — `postprocess_*.py`, `extract_5mer.py`,
`r2d_liftover.py`, `Epinano_DiffErr.R`, the report and plotting helpers — are part
of this benchmark, written to normalise each tool's output into the common
BED-like format.

Every command line that was run against these tools, including the paths of the
tool distributions used, is in
[`../tools/inventory/tables/TI2_command_lines.csv`](../tools/inventory/tables/TI2_command_lines.csv),
with versions in `TI3_software_versions.csv` and model/checkpoint files in
`TI4_model_checkpoints.csv`. See [tool_inventory_notes.md](tool_inventory_notes.md)
for how those were collected and what "static_script" means for their evidence
strength.
