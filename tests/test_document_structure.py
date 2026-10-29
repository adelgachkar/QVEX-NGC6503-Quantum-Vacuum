"""Structural guards for the Markdown record.

These tests do not read the model; they police the documents themselves. Each one
corresponds to a defect that was present in this bundle and that no equation-level
test could see:

* a code fence left open at end of file (``01_Vacuum_Popcorn_Density_Node.md``),
* two documents spliced into one file, each declaring the same identifier
  (``Annex_A_Vacuum_Popcorn_Physical_Constraints.md``),
* an annex naming a baseline module that this bundle does not ship.

They fail loudly rather than skipping, so a document that loses its structure cannot
pass by being unparsable.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Directories that hold third-party, archived or generated material.
EXCLUDED_DIRS = {".venv", "archive", ".pytest_cache", ".obsidian", ".freebuff", ".git", "docs", "outputs"}

DOCUMENT_ID = re.compile(r"\*\*Document ID:\*\*\s*`?([A-Za-z0-9][A-Za-z0-9._\-]*)`?")

ANNEX_A = "nodes_analysis/Annex_A_Vacuum_Popcorn_Physical_Constraints.md"
BUNDLE_NODE = "nodes/active_tensor_foam.py"
PARENT_NODE = "qvex/nodes/popcorn_vacuum_density.py"


def documents():
    """Every Markdown file that belongs to this record."""
    return sorted(
        path for path in ROOT.rglob("*.md")
        if not (set(path.relative_to(ROOT).parts) & EXCLUDED_DIRS)
    )


def read(path):
    return path.read_text(encoding="utf-8")


def fence_lines(text):
    return [i for i, line in enumerate(text.split("\n"), 1) if line.startswith("```")]


def headings_outside_fences(text):
    """H1 headings that are real Markdown, not comments trapped in a code block."""
    found, in_fence = [], False
    for line in text.split("\n"):
        if line.startswith("```"):
            in_fence = not in_fence
            continue
        if not in_fence and line.startswith("# "):
            found.append(line)
    return found


def test_every_document_has_balanced_code_fences():
    """An odd fence count means the rest of the file renders as code."""
    broken = {}
    for path in documents():
        fences = fence_lines(read(path))
        if len(fences) % 2:
            broken[str(path.relative_to(ROOT))] = f"odd fence count {len(fences)} at lines {fences}"
    assert not broken, f"unbalanced code fences: {broken}"


def test_document_identifiers_are_declared_once_and_globally_unique():
    """One identifier per file, and no identifier shared by two files."""
    seen, problems = {}, []
    for path in documents():
        rel = str(path.relative_to(ROOT))
        ids = DOCUMENT_ID.findall(read(path))
        if len(ids) > 1:
            problems.append(f"{rel} declares {len(ids)} identifiers: {ids}")
        for identifier in ids:
            if identifier in seen and seen[identifier] != rel:
                problems.append(f"identifier {identifier} is shared by {seen[identifier]} and {rel}")
            seen[identifier] = rel
    assert not problems, f"identifier problems: {problems}"


def test_governance_documents_have_exactly_one_title():
    """A file that declares an identifier is one document, so it has one real H1.

    Heading counting is fence-aware on purpose: an unclosed fence turns the rest of a
    file into code, which is how a PowerShell comment was once counted as a title.
    """
    problems = {}
    for path in documents():
        text = read(path)
        if not DOCUMENT_ID.search(text):
            continue
        titles = headings_outside_fences(text)
        if len(titles) != 1:
            problems[str(path.relative_to(ROOT))] = titles
    assert not problems, f"documents must carry exactly one H1 outside code fences: {problems}"


def test_annex_a_discloses_that_its_baseline_is_not_shipped_here():
    """The annex governs a parent-engine node; silence about that is a provenance gap."""
    annex = read(ROOT / ANNEX_A)
    bundle_node = ROOT / BUNDLE_NODE

    assert PARENT_NODE in annex, "Annex A no longer names its baseline module"
    assert not (ROOT / PARENT_NODE).exists(), (
        f"{PARENT_NODE} is now shipped in this bundle: the provenance note in Annex A "
        "must be updated, because it states the opposite"
    )
    assert bundle_node.exists(), f"the bundle's executable node {BUNDLE_NODE} is missing"
    assert "is **not shipped**" in annex, (
        "Annex A must state explicitly that its baseline module is not shipped in this bundle"
    )
    assert "ActiveTensorFoamSolver" in annex, (
        "Annex A must name the bundle's actual executable node"
    )
