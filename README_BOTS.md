# Project Layout
Tasks are in `tasks/`. They get imported by the workflows in `workflows/`.
# Command Pallette
* Check syntax: `miniwdl check workflows/<workflow>.wdl`
* Lint: `sprocket lint workflows tasks tests`
* Refresh lint baseline: `sprocket lint --generate-baseline workflows tasks tests`
* Check that inputs, `parameter_meta`, README, and `params/` agree: `python3 scripts/lint_wdl_docs.py`
* Run one task on example inputs: `miniwdl run --as-me tasks/<file>.wdl --task <task> -i params/<inputs>.json`
# Linting
`sprocket.toml` configures lint rules; `sprocket-baseline.toml` records pre-existing
findings so only new ones fail CI. Baseline entries are matched on a hash of the
flagged source, so editing a previously flagged line makes its entry "stale" and
fails the build. Regenerate the baseline and commit it when that happens.
# Project Rules
* Make sure workflow parameters have documentation in the meta section of the workflow.
* Make sure the documentation in the meta section of the workflow agrees with the documentation in the README.
