---
title: "Technical Review and Audit Report: QVEX-NGC6503 Independent Analysis (v2)"
author: "Adel Gachkar"
project: "QVEX-NGC6503-Quantum-Vacuum"
document_type: "Review Report"
version: "2"
supersedes: "QVEX-REV-NGC6503-2026-V1"
date: "2026-10-27"
tags: [QVEX, Verification, CodeReview, NGC6503, FileRanking]
---

# Technical Review and Audit Report: QVEX-NGC6503 Independent Analysis (v2)

**Document ID:** QVEX-REV-NGC6503-2026-V2  
**Supersedes:** [[QVEX_NGC6503_Review_Report_v1|QVEX-REV-NGC6503-2026-V1]] (its stamp described an earlier snapshot)  
**Target:** `docs/NGC6503-Quantum-Vacuum` in `https://github.com/adelgachkar/QVEX_v0.1.0`  
**Local source of truth:** `QVEX-NGC6503-Quantum-Vacuum` (byte-identical copy, verified)  
**Review date:** 2026-10-27  
**Scope:** Re-audit of the corrected package — test execution, structural integrity, cross-document consistency, reproducibility, publication state, and a per-file quality ranking.

---

## 1. Executive Summary & Verification Stamp

| Check | Result |
| :--- | :--- |
| Test command | `PYTHONPATH=. python -m pytest -q` |
| Outcome | **6 passed** · `PYTEST_EXIT = 0` |
| Environment | Windows · Python 3.12.6 · numpy 2.4.0 · pytest 9.1.1 |
| Mutation testing | 3/3 injected defects detected (`exit = 1`); restored file `exit = 0` with identical SHA-256 |
| Code fences | **9/9 balanced** (0 or 2 ` ``` ` lines per document) |
| YAML frontmatter | **9/9 present**, `yaml.safe_load` → `invalid: 0` |
| Wikilinks | **27 total, `broken: 0`** |
| Plot references | single target `outputs/plots/sparc_ngc6503_fit.png` → **RESOLVES**; no stale `_pid_fit` outside this report family |
| Local markdown links | `ok = 1`, `broken = 0` |
| Duplicate filenames | **0** (the two divergent `NGC6503_QVEX_Complete_Analysis.md` copies were merged) |
| Publication | `main = e467944` pushed as `1b8dbc9..e467944`; all 10 files present on `origin/main`; raw files verified HTTP 200 |
| Zenodo metadata | `CITATION.cff` and `.zenodo.json` already carry the supplied identity (Gachkar, Adel · Islamic Azad University, Urmia · ORCID 0009-0006-7713-6004 · adelgachkar@gmail.com) |

**Overall evaluation:** **Approved with reservations.** The numerical pipeline executes, the node/report/test triangle is consistent for the vacuum-density model, and the package is structurally sound and published. Four reservations remain open (§4).

---

## 2. Status of v1 Findings

| v1 finding | Status in this revision | Evidence |
| :--- | :--- | :--- |
| 1. Vacuum density formulation mismatch (node 01 vs test/report) | **Resolved** | Node 01 now carries the canonical equation, byte-identical to `NGC6503_QVEX_Complete_Analysis.md` §1; parameters `γ_vac = 1.0`, `R_cutoff = 50.0` documented in the node |
| 2. Nonlinear foam friction divergence (node 02 vs test) | **Resolved** | Node 02 presents the radial form as canonical and labels the velocity-gradient expression a conceptual precursor |
| 3. Broken image references (`sparc_ngc6503_pid_fit.png`) | **Resolved** | Both references retargeted to the existing `outputs/plots/sparc_ngc6503_fit.png`; `git check-ignore` confirms the committed copy is not ignored |
| 4. Structural artifacts (unclosed fence, zero-byte duplicate) | **Resolved** | Fence closed; the fixed `NGC6503_QVEX_Complete_Analysis.md` (3143 B) replaced the stale root copy and the empty `outputs/plots/` duplicate was removed |
| — stamp items `4 passed` / `Git clean` / zero-byte path | **Superseded** | They described the previous snapshot; v2 replaces them with the measurements above |

---

## 3. Verification Detail

### 3.1 Test suite is bound to Node 01

The vacuum-density test no longer exercises an internal copy of the formula. It parses `$$ … $$` out of `nodes_analysis/01_Vacuum_Popcorn_Density_Node.md`, translates the LaTeX (failing loudly on untranslated commands), parses the documented parameters from the same file, and evaluates the node's own equation against `rho_popcorn()` over 201 radii (`rtol = 1e-12`), plus monotonicity, positivity, and boundary assertions.

### 3.2 Mutation testing (defect injection, all reverted)

| Mutation applied to Node 01 | pytest exit | Interpretation |
| :--- | :---: | :--- |
| Exponent sign flip `^{-γ}` → `^{+γ}` | **1** (2 failed) | numeric comparison **and** report-equation sync both catch it |
| Documented `R_cutoff` 50.0 → 40.0 | **1** (2 failed) | parameter parsing binds the doc to the code defaults |
| Legacy pre-unification formula restored | **1** (1 failed) | the old node/report divergence is detected |
| File restored (SHA-256 identical) | **0** (6 passed) | clean state verified after every mutation |

### 3.3 Structural and link integrity

- **Fences:** every document has an even number of ` ``` ` lines (master note 2 after closing the parameter block; all others 0).
- **Frontmatter:** present and parseable in all 9 markdown documents.
- **Wikilinks:** 27 links, all targets resolve to existing note stems; the master note links all four nodes plus the report, the technical reference and this review family.
- **Plot asset:** the only referenced image resolves relative to each document; the single markdown image embed returns `ok`.

---

## 4. Open Findings (reservations)

| # | Severity | Finding | Suggested action |
| :--- | :--- | :--- | :--- |
| O1 | High | The foam-friction function changes sign at **r ≈ 12.1 kpc** (`F(15) = -0.109`, `F(20) = -0.348`), i.e. it becomes anti-damping over the outer disk, while node 02 claims "stabilization beyond". | Clip/revise the term or document the sign change explicitly; add a regression test |
| O2 | Medium | Tests for nodes 02–04 are still self-referential (formulas re-typed inside the test module). | Apply the node-parsing pattern used for node 01 |
| O3 | Medium | `QVEX_NGC6503_Technical_Reference.md` claims fitting routines are "verified against empirical SPARC datasets", but the package ships no SPARC data and no plot script — the PNG is not reproducible. | Add the data file and the plotting script, or soften the claim |
| O4 | Low | Master note title uses `##` instead of `#`; technical reference embeds personal contact data (email/ORCID/affiliation) — deliberate for citation, worth confirming before wider distribution. | Normalize the heading; confirm the PII decision |

Also noted: the standalone package carries no `.git` (its history lives in `QVEX_v0.1.0`), and `pyproject.toml` restricts `testpaths`, so `docs/NGC6503-Quantum-Vacuum/tests` is not collected by CI unless the path is added.

---

## 5. File Ranking

Criteria (0–10): internal consistency 30% · reproducibility/testability 30% · document structure 20% · quantitative and citation rigor 20%.

| Rank | File | Score | Verdict |
| :---: | :--- | :---: | :--- |
| 1 | `tests/test_ngc6503_vacuum_nodes.py` | **8.0** | Only executable artifact; now parses node 01 from markdown and was mutation-proven. Gap: nodes 02–04 still tested by internal copy |
| 2 | `nodes_analysis/01_Vacuum_Popcorn_Density_Node.md` | **7.5** | Canonical equation, all four parameters documented, frontmatter, linked — previously the weakest note |
| 3 | `NGC6503_QVEX_Complete_Analysis.md` | **7.0** | Single merged copy, §5 completed with a rendered image, linked to all nodes |
| 4 | `nodes_analysis/04_Dipole_Axis_Constraint_Node.md` | **7.0** | Best-structured note; mappings explicitly labelled untested hypotheses |
| 5 | `nodes_analysis/02_Nonlinear_Foam_Friction_Node.md` | **6.5** | Divergence resolved, but the outer-disk sign change (O1) contradicts "stabilization" |
| 6 | `QVEX_NGC6503_Axial_Constraint_Master_Note.md` | **6.5** | Fence closed, plot path fixed, 7 working links — heading level still `##` |
| 7 | `nodes_analysis/03_PID_Controller_Dynamic_Tuning_Node.md` | **6.0** | Gains consistent; stability/phase-margin and the 11.8 kpc "edge" remain unsupported |
| 8 | `QVEX_NGC6503_Technical_Reference.md` | **6.0** | Cleanest prose and correct SPARC citation; over-claim (O3) and embedded PII (O4) |
| 9 | `outputs/plots/sparc_ngc6503_fit.png` | **5.5** | Valid PNG (1350×1200, 96 KB) and correctly referenced, but not reproducible from the package |
| 10 | `review/QVEX_NGC6503_Review_Report_v1.md` | **4.5** | Structurally valid but describes a superseded snapshot; retained for provenance only |

*(This v2 report is not ranked.)*

---

## 6. Limitations of this Review

- `git status` cannot be asserted for the standalone package (no `.git`); publication was verified through `git ls-remote`, `git ls-tree origin/main`, and HTTP 200 fetches of the raw files instead.
- The PNG's scientific content was not cross-checked against SPARC values because neither the data nor the plotting script ships with the package.
- GitHub Actions was not executed; only local collection (`65 tests collected`, `COLLECT_EXIT = 0`) was observed.
- Zenodo publication and DOI issuance are platform actions and were not performed.

## 7. Recommendations

1. Resolve O1 (friction sign) before any physical interpretation of the outer-disk fit.
2. Extend the node-parsing test pattern to nodes 02–04 so every documented equation is machine-checked.
3. Ship the SPARC input and the plotting script to make `sparc_ngc6503_fit.png` reproducible (O3).
4. Add `docs/NGC6503-Quantum-Vacuum/tests` to `testpaths` so CI enforces the node/report/test triangle.
