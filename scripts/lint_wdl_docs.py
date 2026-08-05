#!/usr/bin/env python3
"""
Check that the workflows and their documentation agree.

README sections are matched to workflows by looking for a link to the workflow
file in a section, and only parameter names are compared; the descriptions may
differ.

Run from the root of the repository, or pass the root as the only argument.
"""

import os
import re
import sys

import WDL

# Where each kind of file lives, relative to the repository root.
WORKFLOW_DIR = "workflows"
# Subworkflows that exist only to be called by other workflows. They are held to
# the parameter_meta rules like anything else, but they are plumbing, so the
# README describes them in a sentence instead of documenting their parameters.
INTERNAL_DIR = os.path.join("workflows", "internal")
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


def find_wdls(directory):
    """
    Return the paths of every .wdl file at or below a directory, sorted. This
    picks up workflows/internal, where the subworkflows that exist only to be
    called by other workflows live.
    """

    found = []
    for parent, _, names in os.walk(directory):
        found += [os.path.join(parent, n) for n in names if n.endswith(".wdl")]
    return sorted(found)


def load_documents(directory):
    """
    Parse every .wdl file at or below a directory. Returns a list of
    (path, document).

    Files that don't parse are skipped: making sense of them is ``miniwdl
    check``'s job, and it runs alongside this.
    """

    documents = []
    for path in find_wdls(directory):
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
    Find the README section that documents a workflow, by looking for the
    workflow's path anywhere in the section. That covers a relative link, an
    absolute github.com link, and a `miniwdl run` example alike. Returns
    (heading, body_lines), or None if no section mentions the workflow.
    """

    for heading, body in sections:
        for line in body:
            if workflow_path in line:
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
    Check a workflow's README section against its inputs.

    A workflow that users are meant to run has to have a section, and that
    section has to list every parameter. Internal subworkflows are skipped: the
    README mentions them, but making callers read a parameter list for plumbing
    they never invoke would be noise.
    """

    if path.startswith(INTERNAL_DIR + os.sep):
        return

    section = readme_section_for(sections, path)
    if section is None:
        problems.report(path, "has no section in {} that mentions it".format(README))
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


def main(args):
    root = args[0] if args else "."
    os.chdir(root)

    problems = Problems()
    exempt = read_exemptions(EXEMPTIONS)

    workflow_documents = load_documents(WORKFLOW_DIR)
    sections = read_readme_sections(README)

    for path, document in workflow_documents:
        if document.workflow is None or path in exempt:
            continue
        check_parameter_meta(problems, path, document.workflow)
        check_readme(problems, path, document.workflow, sections)

    if problems.count:
        print("\n{} documentation problem(s) found".format(problems.count), file=sys.stderr)
        return 1
    print("Workflows and their documentation agree")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
