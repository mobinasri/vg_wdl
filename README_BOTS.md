# Project Layout
Tasks are in `tasks/`. They get imported by the workflows in `workflows/`.
# Command Pallette
* Check syntax: `miniwdl check workflows/<workflow>.wdl`
* Lint: `sprocket lint workflows tasks tests`
* Check that inputs, `parameter_meta`, README, and `params/` agree: `python3 scripts/lint_wdl_docs.py`
* Run one task on example inputs: `miniwdl run --as-me tasks/<file>.wdl --task <task> -i params/<inputs>.json`
# Linting
`sprocket.toml` configures lint rules. Notes are non-fatal; only warnings fail
CI. Rules with widespread pre-existing findings are excepted wholesale in
`sprocket.toml` rather than tracked in a baseline file, so any lint failure is
either a new finding in code you touched or an intentional `#@ except`
directive to add.
# Project Rules
* Make sure workflow parameters have documentation in the meta section of the workflow.
* Make sure the documentation in the meta section of the workflow agrees with the documentation in the README.
