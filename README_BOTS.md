# Project Layout
Tasks are in `tasks/`. They get imported by the workflows in `workflows/`.
# Command Pallette
* Check syntax: `miniwdl check workflows/<workflow>.wdl`
* Lint: `sprocket lint workflows tasks tests`
* Refresh lint baseline: `sprocket lint --generate-baseline workflows tasks tests`
# Linting
`sprocket.toml` configures lint rules; `sprocket-baseline.toml` records pre-existing
findings so only new ones fail CI. Baseline entries are matched on a hash of the
flagged source, so editing a previously flagged line makes its entry "stale" and
fails the build. Regenerate the baseline and commit it when that happens.
