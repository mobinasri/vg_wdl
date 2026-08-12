# Project Layout
Tasks are in `tasks/`. They get imported by the workflows in `workflows/`.
# Command Pallette
* Check syntax: `miniwdl check workflows/<workflow>.wdl`
* Check that inputs, `parameter_meta`, README, and `params/` agree: `python3 scripts/lint_wdl_docs.py`
* Run one task on example inputs: `miniwdl run --as-me tasks/<file>.wdl --task <task> -i params/<inputs>.json`
# Project Rules
* Make sure workflow parameters have documentation in the meta section of the workflow.
* Make sure the documentation in the meta section of the workflow agrees with the documentation in the README.
