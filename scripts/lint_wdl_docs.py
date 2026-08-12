#!/usr/bin/env python3
"""
Check that the workflows and their documentation agree.

README sections are matched to workflows by looking for a reference to the
workflow file in a section, and only parameter names are compared; the
descriptions may differ.

Run from the root of the repository, or pass the root as the only argument.
"""

import os
import re
import sys
from pathlib import Path
from typing import Iterator

import WDL

# This is the README file to check workflows against
README = "README.md"
# This directory has all the workflows 
WORKFLOW_DIR = Path("workflows")
# Subworkflows that exist only to be called by other workflows.
# We skip these.
INTERNAL_DIR = WORKFLOW_DIR / "internal"
# We also skip the workflows listed in this file.
EXEMPTIONS = Path("scripts") / "doc_lint_exemptions.txt"

# Regex to match a parameter entry in th README: "- *NAME*: description".
README_PARAMETER_RE = re.compile(r"^\s*[-*]\s+\*([A-Za-z_][A-Za-z_0-9]*)\*\s*:")
# Regex to match a Markdown section heading line
README_HEADING_RE = re.compile(r"^(#+)\s+(.*\S)\s*$")


class Report:
    """
    Class for counting and reporting lint problems.
    """

    def __init__(self) -> None:
        self.problems = 0

    def problem(self, where: str, message: str) -> None:
        """
        Report a lint problem at a particular place with a particular message.
        """
        print("{}: {}".format(where, message), file=sys.stderr)
        self.problems += 1

class Readme:
    """
    Represents the README and allows pulling documented workflow inputs from it.
    """

    def __init__(self, path: Path) -> None:
        """
        Parse the README.
        """

        # This holds README sections as lists of lines.
        self.sections = []
        for line in open(path):
            heading = README_HEADING_RE.match(line)
            if heading:
                # Start a new section
                self.sections.append([line])
            elif self.sections:
                # At least one section exists, so one is in progress.
                self.sections[-1].append(line)

    def get_section(self, workflow_filename: str) -> list[str] | None:
        """
        Find the README section that documents a workflow.

        Looks for the workflow's path anywhere in the section.
        """

        for lines in self.sections:
            # TODO: It's not super fast to scan each section all the time. We
            # should maybe make a rule about the headings.
            for line in lines:
                if workflow_filename in line:
                    return lines
        return None

    def get_inputs(self, workflow_filename: str) -> set[str] | None:
        """
        Get the inputs for a workflow, as listed in the README.
        """

        section = self.get_section(workflow_filename)

        if section is None:
            # Not found
            return None

        parameters = set()

        for line in section:
            match = README_PARAMETER_RE.match(line)
            if match:
                parameters.add(match.group(1))

        return parameters


def read_set(path: Path) -> set[str]:
    """
    Read a set of strings from a file.

    Ignores # comments.
    """
    
    return {l.split("#", 1)[0].strip() for l in open(path).readlines()} - {""}


def scan_workflows(report: Report) -> Iterator[WDL.Tree.Workflow]:
    """
    Get an iterator over WDL workflows we want to process.

    Remember that WDL workflows have a .pos.uri for where they came from.
    """

    skip_set = read_set(EXEMPTIONS)

    for path in WORKFLOW_DIR.rglob("*.wdl"):
        if str(path) in skip_set:
            # Ignore listed ignored files
            continue
        if INTERNAL_DIR in path.parents:
            # Ignore all the internal workflows
            continue
        try:
            document = WDL.load(str(path))
        except Exception as e:
            report.problem(path, f"does not parse: {e}")
            continue

        if document.workflow:
            yield document.workflow


def check_workflow(report: Report, readme: Readme, workflow: WDL.Tree.Workflow) -> None:
    """
    Find lint errors in a workflow.
    """
    
    # Get all the workflow input names
    workflow_inputs = {d.name for d in (workflow.inputs or [])}

    # Make sure the input names are 1 to 1 with the README documentation entries
    readme_inputs = readme.get_inputs(workflow.pos.uri)
    if readme_inputs is None:
        report.problem(workflow.pos.uri, f"has no section in README")
    else:
        for name in sorted(workflow_inputs - readme_inputs):
            report.problem(workflow.pos.uri, f"input {name} is not listed in the README")
        for name in sorted(readme_inputs - workflow_inputs):
            report.problem(workflow.pos.uri, f"README lists {name}, which is not an input")


def main(args: list[str]) -> int:
    root = args[0] if args else "."
    os.chdir(root)

    report = Report()

    readme = Readme(README)

    for workflow in scan_workflows(report):
        check_workflow(report, readme, workflow)

    if report.problems:
        print("\n{} problem(s) found".format(report.problems), file=sys.stderr)
        return 1
    print("Workflow documentation OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
