#!/usr/bin/env python3
"""
Check that the workflows, their documentation, and the example inputs agree.

Three things drift apart as workflows are edited, and none of them are caught by
``miniwdl check``: the ``parameter_meta`` section stops listing every input, the
README stops listing every parameter, and the example inputs in ``params/`` keep
referring to inputs that have been renamed or removed. This finds all three.

The README is matched to a workflow by looking for a link to the workflow file in
a section, and only parameter *names* are compared, so the README is free to
reword a description without failing here.

Run from the root of the repository, or pass the root as the only argument.
"""

import json
import os
import re
import sys

import WDL

# Where each kind of file lives, relative to the repository root.
WORKFLOW_DIR = "workflows"
TASK_DIR = "tasks"
PARAMS_DIR = "params"
README = "README.md"
# Files that are known not to pass and that nobody is expected to fix. See the
# file itself for what belongs in it.
EXEMPTIONS = os.path.join("scripts", "doc_lint_exemptions.txt")

# A parameter as the README lists it: "- *NAME*: description".
README_PARAMETER_RE = re.compile(r"^\s*[-*]\s+\*([A-Za-z_][A-Za-z_0-9]*)\*\s*:")
# A Markdown section heading, at any depth.
README_HEADING_RE = re.compile(r"^(#+)\s+(.*\S)\s*$")


class Problems:
    """
    Collects problems as they are found, so that one run reports all of them
    instead of stopping at the first.
    """

    def __init__(self):
        self.count = 0

    def report(self, where, message):
        print("{}: {}".format(where, message), file=sys.stderr)
        self.count += 1


def read_exemptions(path):
    """
    Read the set of exempt file paths. Blank lines and #-comments are ignored.
    """

    exempt = set()
    if not os.path.exists(path):
        return exempt
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                exempt.add(line)
    return exempt


def load_documents(directory):
    """
    Parse every .wdl file in a directory. Returns a list of (path, document).

    Files that don't parse are skipped: making sense of them is ``miniwdl
    check``'s job, and it runs alongside this.
    """

    documents = []
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".wdl"):
            continue
        path = os.path.join(directory, name)
        try:
            documents.append((path, WDL.load(path)))
        except Exception as e:
            # Only the first line: parse errors from the WDL grammar list every
            # token that would have been accepted, which is pages of noise.
            summary = str(e).splitlines()[0] if str(e) else type(e).__name__
            print("{}: skipped, does not parse: {}".format(path, summary), file=sys.stderr)
    return documents


def declared_inputs(workflow):
    """
    Return the names of a workflow's own inputs.

    Only the input section counts. Inputs that a called subworkflow leaves
    unbound are also settable from outside, but they belong to the subworkflow's
    documentation, not to this one's.
    """

    return [declaration.name for declaration in (workflow.inputs or [])]


def required_inputs(workflow):
    """
    Return the names of the inputs that a caller has to supply: the ones with no
    default and no way to be left out.
    """

    required = []
    for declaration in workflow.inputs or []:
        if declaration.expr is None and not declaration.type.optional:
            required.append(declaration.name)
    return required


def read_readme_sections(path):
    """
    Return a list of (heading, body_lines) for the README, one per heading, in
    document order. Text before the first heading is dropped.
    """

    sections = []
    with open(path) as f:
        for line in f:
            heading = README_HEADING_RE.match(line)
            if heading:
                sections.append((heading.group(2), []))
            elif sections:
                sections[-1][1].append(line)
    return sections


def readme_section_for(sections, workflow_path):
    """
    Find the README section that documents a workflow, by looking for a link to
    the workflow's file. Returns (heading, body_lines), or None if no section
    claims the workflow.
    """

    needle = "(" + workflow_path + ")"
    for heading, body in sections:
        for line in body:
            if needle in line:
                return (heading, body)
    return None


def check_parameter_meta(problems, path, workflow):
    """
    Every input needs a parameter_meta entry, and every parameter_meta entry
    needs an input, so that the documentation can't quietly describe a parameter
    that no longer exists.
    """

    inputs = set(declared_inputs(workflow))
    documented = set(workflow.parameter_meta.keys())

    for name in sorted(inputs - documented):
        problems.report(path, "input {} has no parameter_meta entry".format(name))
    for name in sorted(documented - inputs):
        problems.report(path, "parameter_meta documents {}, which is not an input".format(name))


def check_readme(problems, path, workflow, sections):
    """
    Check the README section for a workflow, if it has one, against the
    workflow's inputs. Workflows with no section are left alone; a lot of the
    older ones have never had one.
    """

    section = readme_section_for(sections, path)
    if section is None:
        return
    heading, body = section

    inputs = set(declared_inputs(workflow))
    documented = set()
    for line in body:
        match = README_PARAMETER_RE.match(line)
        if match:
            documented.add(match.group(1))

    where = "{} ({} section of {})".format(path, heading, README)
    for name in sorted(inputs - documented):
        problems.report(where, "input {} is not listed in the README".format(name))
    for name in sorted(documented - inputs):
        problems.report(where, "README lists {}, which is not an input".format(name))


def check_params_file(problems, path, callables):
    """
    Check one example inputs file: every key has to name a real input of the
    workflow or task it is namespaced under, and every required input of that
    workflow or task has to be set.

    ``callables`` maps a workflow or task name to its object.
    """

    with open(path) as f:
        try:
            values = json.load(f)
        except ValueError as e:
            problems.report(path, "is not valid JSON: {}".format(e))
            return

    for key in values:
        if "." not in key:
            problems.report(path, "key {} is not of the form NAME.INPUT".format(key))
            continue
        target_name, input_name = key.split(".", 1)
        target = callables.get(target_name)
        if target is None:
            problems.report(path, "key {} names no workflow or task".format(key))
            continue
        if "." in input_name:
            # An override aimed at an input of a call inside the workflow, like
            # WORKFLOW.someCall.some_input. Which calls exist is the workflow's
            # business and can change, so we only check the workflow name.
            continue
        if input_name not in declared_inputs(target):
            problems.report(path, "{} is not an input of {}".format(input_name, target_name))

    # Only complain about missing inputs for things this file actually sets up,
    # so that a file for one workflow says nothing about any other.
    for target_name in sorted({key.split(".", 1)[0] for key in values if "." in key}):
        target = callables.get(target_name)
        if target is None:
            continue
        for name in required_inputs(target):
            if "{}.{}".format(target_name, name) not in values:
                problems.report(path, "required input {}.{} is not set".format(target_name, name))


def main(args):
    root = args[0] if args else "."
    os.chdir(root)

    problems = Problems()
    exempt = read_exemptions(EXEMPTIONS)

    workflow_documents = load_documents(WORKFLOW_DIR)
    sections = read_readme_sections(README)

    # Everything a params file could be for: the workflows, plus the tasks, so
    # that a task can be given example inputs and run on its own.
    callables = {}
    for _, document in workflow_documents:
        if document.workflow is not None:
            callables[document.workflow.name] = document.workflow
    for _, document in load_documents(TASK_DIR):
        for task in document.tasks:
            callables[task.name] = task

    for path, document in workflow_documents:
        if document.workflow is None or path in exempt:
            continue
        check_parameter_meta(problems, path, document.workflow)
        check_readme(problems, path, document.workflow, sections)

    for name in sorted(os.listdir(PARAMS_DIR)):
        path = os.path.join(PARAMS_DIR, name)
        if name.endswith(".json") and path not in exempt:
            check_params_file(problems, path, callables)

    if problems.count:
        print("\n{} documentation problem(s) found".format(problems.count), file=sys.stderr)
        return 1
    print("Workflow documentation and example inputs agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
