# Test Verification Summary

- Canonical suite: `tests/` (see `pytest.ini` → `testpaths = tests`).
- `python -m pytest -q`: **32 passed** (4 modules + 1 helper, exit 0).
- `python -m pytest -q .` (whole tree): **32 passed**, exit 0 — no collection errors.
- Environment: Python 3.12.6 · numpy 2.4.0 · pytest 9.1.1 (`-W error::UserWarning` clean).

## Document ↔ code binding (Nodes 01, 02 and 05)

`tests/node_parser.py` extracts the display equations from the node Markdown, translates
the LaTeX subset (`\frac`, `^{...}`, `\text{...}`, `\exp`, `\ln`, `\max`, implicit
multiplication) into Python and evaluates it with the implementation's own parameters.
Untranslatable LaTeX raises instead of being skipped.

| Test | Binds |
| :--- | :--- |
| `test_node_01_density_matches_the_documented_equation` | Node 01 formula form vs `vacuum_density()`, rtol 1e-12 over 201 radii |
| `test_node_01_documented_parameters_match_the_implementation` | documented `lambda_coh`, `gamma_vac`, `R_cutoff`, `rho_0` vs code defaults |
| `test_node_02_raw_friction_matches_the_documented_equation` | Node 02 formula vs `foam_friction_raw()`, rtol 1e-12 |
| `test_node_02_clamp_is_removed_and_the_production_law_never_reverses` | `foam_friction` stays deleted; the exchange force is strictly damping for `r > 0` |
| `test_node_02_sign_change_is_located_and_documented` | the documented expression still reverses at 12.0 < r < 12.2 kpc |
| `test_node_05_exchange_coefficient_matches_the_documented_equation` | Node 05 formula vs `exchange_coefficient()`, rtol 1e-12 |
| `test_node_05_documented_parameters_match_the_implementation` | documented `Gamma_0`, `r_core`, `lambda_coh` vs code defaults |
| `test_node_05_force_sign_is_set_by_relative_motion_only` | no radius can reverse the force at fixed relative velocity |
| `test_node_05_matter_and_vacuum_terms_sum_to_zero` | the exchange pair is conservative exactly |

The document is the source of truth: changing either side alone turns the suite red.

Mutation evidence (each reverted, file SHA-256 verified identical afterwards):

| Injected defect | Result |
| :--- | :--- |
| code `gamma_vac` 1.0 → 1.2 | exit 1 (parameter binding caught it) |
| code `exp(-r/R_cutoff)` → `exp(-r/2R_cutoff)` | exit 1 (form binding caught it) |
| document Node 02 amplitude 1.1 → 1.2 | exit 1 (formula binding caught it) |

## Model changes (2026-10-29)

- `vacuum_density` is now Node 01's canonical profile
  `rho_0 (1 + r/lambda_coh)^(-gamma_vac) exp(-r/R_cutoff)`, amplitude `7.02e-30 g/cm^3`
  as documented (previously `rho_0 = 0.045` unitless with a `(r/r_core)^1.2` law).
- `foam_friction_raw` is Node 02's documented expression verbatim, sign included.
- The `max(0, raw)` clamp has been **removed**, not kept: `foam_friction` no longer exists
  and `compute()` returns `exchange_coefficient` and `exchange_friction` instead.
- Known limitation, documented in code: `effective_enclosed_mass` still mixes
  `g/cm^3` with a kpc/Msun `G` and is not bound by any test.

## Deduplication (2026-10-29)

Five files with duplicate or dead content were retired to `archive/2026-10-29-dedup/`
(excluded from collection via `norecursedirs`):

| Retired file | Reason |
| :--- | :--- |
| `root__test_active_tensor_foam.py` | duplicate of `tests/test_active_tensor_foam.py` |
| `root__test_ngc6503_vacuum_nodes.py` | duplicate of the `nodes_analysis/` copy |
| `nodes_analysis__test_ngc6503_vacuum_nodes.py` | same module under a second basename; merged into `tests/` |
| `nodes_analysis__test_active_tensor_foam.py` | unimportable (referenced `FoamParameters` and methods that do not exist) |
| `tests__active_tensor_foam.py` | dead second implementation of `ActiveTensorFoamNode`; nothing imported it |

`test_only_one_test_module_per_basename_outside_the_archive` fails if duplicate
`test_*.py` basenames reappear (verified by injecting one).

The `.venv` in this bundle ships with `pip` only (no numpy/pytest) and the project has no
`requirements.txt` or `pyproject.toml`; the runs above used the system interpreter.

## Exchange current replaces the clamped friction (2026-10-29)

Node 02 changes sign near 12.11 kpc because of its penalty term, so any production law
built on it needs a clamp. The clamp hid the defect instead of removing it. The drag law is
now the Node 05 causal exchange current, `J^mu = Gamma (u_v^mu - u_m^mu)`, with
`F = -Gamma * du`. The coefficient is the Node 02 activation profile without the penalty:
non-negative by construction, so the sign of the force is carried by the relative velocity
alone and no guard is needed.

| r (kpc) | `Gamma` | `F` at du = 1 | `u_m·J` | Node 02 raw |
| ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4.50 | 0.2384 | -0.2384 | -0.2384 | +0.2384 |
| 11.00 | 0.8165 | -0.8165 | -0.8165 | +0.0365 |
| 12.11 | 0.9129 | -0.9129 | -0.9129 | -0.0003 |
| 20.00 | 1.5120 | -1.5120 | -1.5120 | -0.3480 |
| 50.00 | 2.7010 | -2.7010 | -2.7010 | -2.7590 |

The documented expression still reverses; it is kept for comparison, not for use. The
production law has one sign at every radius, and `u_m·J <= 0` everywhere, so the exchange
can only remove energy from the matter sector.

Mutation evidence (each reverted; every file byte-identical afterwards):

| Injected defect | Result |
| :--- | :--- |
| penalty term put back inside `Gamma` | exit 1 (single-sign test caught it) |
| force sign flipped in `exchange_friction` | exit 1 (single-sign test caught it) |
| documented `Gamma_0` 1.1 → 1.2 | exit 1 (parameter binding caught it) |
| clamp `foam_friction` reintroduced | exit 1 (removal test caught it) |
| radius validation deleted | exit 1 (negative-radius test caught it) |

Consequences recorded in the documents: Node 02 now states that its penalty does **not**
ensure non-negativity, and `NGC6503_QVEX_Complete_Analysis.md` §2 no longer calls that
expression the canonical production model.

## Structural repair (2026-10-29)

Four document-level defects were repaired. None of them was visible to an equation-level
test, which is why none had been caught.

| Defect | Before | Repair |
| :--- | :--- | :--- |
| `Annex_A_Vacuum_Popcorn_Physical_Constraints.md` carried two documents | The second annex had been appended inside an A.5 bullet (line 85) of the first, so one file declared `QVEX-ANNEX-A-POP-01` twice with two baselines and two epistemic labels; the `powershell` fence opened at line 159 was never closed | Split at the splice point into Annex A (`QVEX-ANNEX-A-POP-01`) and `nodes_analysis/Annex_A2_Vacuum_Popcorn_Boundary_Conditions.md` (`QVEX-ANNEX-A-BC-01`); the fence was closed; the interrupted A.5 list is restored to six contiguous bullets |
| `01_Vacuum_Popcorn_Density_Node.md` ended inside a fence | 21 lines, one fence: the parameter block never closed | Closing fence added; 22 lines, two fences |
| Annex A named a baseline this bundle does not ship | `qvex.nodes.popcorn_vacuum_density.PopcornVacuumDensityNode` lives in QVEX v0.1.0, yet A.2 and A.4.3 report its `lambda_coh_kpc` and `vacuum_energy_density` as inactive | Added a provenance note in A.1: the module is not shipped here, the bundle node is `ActiveTensorFoamSolver`, and `lambda_coh`/`rho_0` are live parameters in this bundle |
| `ngc6503_metrics.json` claimed verification it does not have | Status read `Verified against SPARC Lelli et al. (2016)` while README states there is no density/friction-to-velocity mapping | Fixed at the generator (`reproduce_ngc6503.py`), not only in the output; the guard line now states that the clamp is numerical |

Text preservation: every technical line of the original Annex A appears verbatim in one of
the two files. The only baseline lines that no longer match are the splice line and the
duplicated identifier — exactly the two declared edits.

Byte-level checks after the repair:

- `ngc6503_residuals.csv` is byte-identical to the pre-repair baseline
  (sha256 `303ac3eb…`), so no number moved: `chi2_total = 58346.27`,
  `chi2_reduced = 1882.14`, `dof = 31` are unchanged.
- `ngc6503_metrics.json` differs from the baseline only in the two status strings.
- `reproduce_ngc6503.py` now writes both artifacts with an explicit LF newline, so the
  outputs no longer depend on the host line-ending convention.

## Structural guards

`tests/test_document_structure.py` binds those defects to the suite:

| Test | Catches |
| :--- | :--- |
| `test_every_document_has_balanced_code_fences` | an unclosed fence in any collected Markdown file |
| `test_document_identifiers_are_declared_once_and_globally_unique` | more than one identifier in a file, or one identifier shared by two files |
| `test_governance_documents_have_exactly_one_title` | a spliced document, counted with fence awareness so a shell comment is not read as a title |
| `test_annex_a_discloses_that_its_baseline_is_not_shipped_here` | an annex that governs a module absent from the bundle without saying so |

Mutation evidence (each reverted; file hash compared afterwards):

| Injected defect | Result |
| :--- | :--- |
| stray fence appended to Node 02 | exit 1 (fence parity failed) |
| `QVEX-ANNEX-A-POP-01` re-inserted into Annex A2 | exit 1 (identifier uniqueness failed) |
| second H1 prepended to Annex A | exit 1 (single-title test failed) |
| provenance phrase `is **not shipped**` removed | exit 1 (disclosure test failed) |

One restore during that run re-encoded Node 02 with LF endings. It was restored from the
untouched baseline copy byte-for-byte (sha256 `8c633774…`); the file is CRLF again and
matches the baseline exactly.
