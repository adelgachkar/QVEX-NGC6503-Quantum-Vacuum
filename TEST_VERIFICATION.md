# Test Verification Summary

- Canonical suite: `tests/` (see `pytest.ini` → `testpaths = tests`).
- `python -m pytest -q`: **88 passed** (8 modules + 2 helpers, exit 0) — 2026-10-30.
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

## One density, one unit (2026-10-30)

The bundle carried the same vacuum amplitude in three unit domains and converted between
them by copying the number:

| Where | Stated value | Unit | Quantity in kg/m^3 |
| :--- | ---: | :--- | ---: |
| Node 01 | 7.02 × 10^-30 | g/cm^3 | 7.02e-27 |
| Annex A2 | ≤ 1e-25 | kg/m^3 | bound only |
| `effective_enclosed_mass` | 7.02e-30 | passed raw into a kpc/Msun solver | wrong by 1.5e31 |

Fix, in three parts:

1. `nodes/active_tensor_foam.py` publishes `CANONICAL_DENSITY_UNIT = "kg/m^3"`, the exact
   SI/IAU constants behind the conversions, `convert_density(value, from, to)`,
   `to_kg_m3`/`from_kg_m3`, and `rho_0_kg_m3` / `rho_0_msun_kpc3`. Unknown units raise
   instead of being guessed. `vacuum_density(r, unit=None)` converts only when asked.
2. `effective_enclosed_mass` converts to `Msun/kpc^3` before meeting `G` in kpc/Msun: the
   halo mass inside 1 kpc is now 2239 Msun instead of 1.5e-28. This closes the limitation
   that the previous revision of this file recorded as "not bound by any test".
3. Node 01 and Annex A2 each carry the same canonical-unit table, and Annex A2 now uses
   Node 01's symbols (`lambda_coh`, `R_cutoff`, `gamma_vac`); the withdrawn `r_s`, `r_c`,
   `alpha` are mapped one-to-one and a test fails if they return to the profile equation.

`tests/test_density_units.py` (9 tests) parses the unit table out of *both* documents,
converts each column to kg/m^3 and compares. The comparison passes `atol=0.0` on purpose:
the amplitudes are ~1e-27, far below NumPy's default `atol=1e-8`, so a default comparison
called a factor-1000 error "close". That was found while proving the test could fail — the
first mutation run passed, which is exactly the failure mode a red test is supposed to
prevent.

Mutation evidence (each reverted; both files byte-identical afterwards):

| Injected defect | Result |
| :--- | :--- |
| Annex A2 `g/cm^3` column overwritten with the `kg/m^3` number | exit 1 (3 tests: cross-document, internal consistency, bound) |
| code `rho_0` 7.02e-30 → 7.02e-27 | exit 1 (document binding + enclosed mass) |

## Non-equilibrium front stress (2026-10-30)

Nodes 01/02/05 are equilibrium laws: `Gamma(r)` is a function of radius, so it cannot
distinguish a front that is arriving from one that has passed, and it has no memory.

`nodes/non_equilibrium_stress.py` implements `Pi^{mu nu} = Delta chi^{mu nu}` in an
orthonormal tetrad, so `u.Pi = 0`, `tr Pi = 0` and symmetry are consequences of the basis
rather than conditions to check after the fact. Israel-Stewart relaxation gives
`Pi(t) = -2 eta sigma (1 - exp(-t/tau))`, i.e. at a front crossing
`Pi = -2 eta sigma (1 - exp(-De))` with `De = tau_Pi v / w`.

One shell `v(r) = (v0/2)[tanh((r-r_in)/w) - tanh((r-r_out)/w)]` produces both fronts from
one field: `sigma_leading = -v0/(3w)` (compression) and `sigma_trailing = +v0/(3w)`
(rarefaction), so the two fronts carry opposite stress. For `tau_Pi = 20 Myr`,
`w = 0.5 kpc`, `v0 = 0.1 kpc/Myr` the delivered stress is 22.1% of the equilibrium value —
a discrepancy no radius-only law can express.

The entropy sign is stated and tested: `Pi : sigma < 0` for a dissipative front, and the
second law is `-(1/2) Pi : sigma = eta sigma : sigma >= 0`. `entropy_production()` returns
the latter.

Causality is enforced, not assumed: `v_c = sqrt(eta / (tau_Pi (rho+p))) <= 1`, so `tau -> 0`
(the instantaneous parabolic law that Nodes 01/02/05 use implicitly) is rejected.

Mutation evidence (each reverted; file hash compared afterwards):

| Injected defect | Result |
| :--- | :--- |
| entropy production sign flipped (`-0.5` → `+0.5`) | exit 1 (3 tests) |
| trailing-front shear sign flipped | exit 1 (opposite-front and dissipation tests) |
| causality guard removed from `validate_transport` | exit 1 (causality test) |

`nodes_analysis/06_Non_Equilibrium_Front_Stress_Node.md` records the algebra, the table of
the three Deborah-number cases, and an explicit list of what is *not* implemented: no
equations of state, no metric or field-equation solve, no angular-momentum sector, no
observational binding, and `T_m` unspecified.

## The exchange four-vector (2026-10-30)

Node 05 states `J^mu = Gamma (u_v - u_m)` and then keeps `F = -Gamma du` and
`u_m.J = -Gamma (du)^2`: a radial, scalar pair. The azimuthal component — the one that
carries torque — was never computed, and Node 06 recorded that as a limitation.

`nodes/exchange_current.py` keeps all four components:
`J^mu = Q u^mu + f^mu` with `u.f = 0`, `Q = -u.J`, `f = h^alpha{}_nu J^nu`.
`ExchangeCurrent` rejects a force with a component along `u`, and `from_vector` splits any
four-vector exactly, so the split is a definition rather than a convention applied later.

Contracting `nabla.T_vac = -J` with `u` gives, term by term, an identity rather than a
model statement:

```
Q = rho_dot + (rho + p_t) Theta + Delta (s^mu s^nu d_mu u_nu) + Pi : sigma
```

with `rho_dot = u.d rho` the convective derivative. Three details are load-bearing: the
pressure in the `p dV` term is `p_t` (not the isotropic mean), the `Delta` term is the
stretching of the anisotropy axis along itself, and `Pi : sigma` survives where
`Delta = 0` because the dissipative stress is its own source.

The projected matter equation `(rho_m + p_m) a^alpha + h^{alpha mu} d_mu p_m = f^alpha` is
implemented as `euler_acceleration()`, and `r f^phi` as `torque_density()`.

Numerical binding: the identity is compared with an independent finite-difference
divergence of `vacuum_stress` on a `(t, r)` grid with a moving front and an azimuthal flow.
Relative residual falls with refinement (`1.5e-2`, `4.8e-3`, `2.7e-3` at 41, 81, 161
points) while the signal stays at `0.16`, so what remains is discretization.

One defect was found by writing the test: the first version projected with the all-up
tensor `h^{mu nu}` applied to a contravariant vector, which does not vanish for `v = u`.
The mixed projector `h^mu{}_nu = delta + u^mu u_nu` is the correct one, and the two are now
separate functions with separate tests. A rest-frame-only implementation hides this error,
which is why it survived until a boosted frame was tested.

Mutation evidence (each reverted afterwards):

| Injected defect | Result |
| :--- | :--- |
| `Pi` dropped from the stress tensor | exit 1 (identity test) |
| `p_t` replaced by the isotropic mean | exit 1 (identity test) |
| projection built from the all-up tensor | exit 1 (6 tests) |
| torque sign reversed | exit 1 (2 tests) |
| `d_t` dropped from the divergence | exit 1 (identity test) |
| `s = e_r` handed to a moving flow | exit 1 (11 tests) |

`nodes_analysis/07_Exchange_Four_Vector_Conservation_Node.md` records the algebra, the
frame and index conventions, the torque table, and the boundary: no equations of state, no
metric, no field-equation solve, no matrix-sector closure, no observational binding.

## The `phi-r` shear channel, and what it can support (2026-10-30)

A recurring claim is that a `Pi^{phi r}` component of the vacuum stress transfers angular
momentum to the gas and thereby holds a flat rotation curve together without dark matter.
The claim decomposes into three assertions with different verdicts, and
`nodes/angular_momentum.py` makes each of them checkable.

| Assertion | Verdict |
| :--- | :--- |
| `Pi^{phi r}` can transport angular momentum | true, but **not** with `Pi = Delta chi`, which has `chi^{r phi} = 0` identically |
| That transport can hold the curve up | false: the `phi` equation holds `T^r{}_phi`, the radial balance holds `p_r`; a torque is not centripetal support |
| So no energy source is needed | false: the torque costs `eta v^2 / r^2` and the same stress dissipates `eta v^2 / (2 r^2)` |

The shear is derived rather than quoted: `sigma_{r phi} = (1/2) r d_r(gamma Omega)` follows
from the coordinate-basis covariant derivative, including the `2/r` connection term that
absence would triple the answer. The suite recomputes it by finite differences
(`rtol 1e-3` at 2000 points) and checks the three standard cases: zero for solid-body
rotation, `-v/(2r)` for a flat curve, `-3 Omega/4` for Keplerian.

Numbers at `r = 10` kpc, `v = 0.1175 kpc/Myr`, `rho_v = 103.7 Msun/kpc^3` and the causality
ceiling `eta = tau rho_v c^2 = 195`: torque density `0.2291`, `P_torque = 0.0269`,
dissipation `0.0135`, total `0.0404 Msun/(kpc Myr^3)`. Divided by `rho_v c^2 = 9.75`, that is
`4.14e-3` per Myr, i.e. an e-folding of 241 Myr against an orbital period of 535 Myr — the
vacuum would be spent within one rotation. The shear stress itself is `2.29`, which is
`6.6e-3` of the gas ram pressure for `rho_gas = 2.5e4 Msun/kpc^3` and `1.7e-3` for `1e5`.
Reaching 10% of that ram pressure would need `eta = 1.2e4`, sixty times the causality
ceiling.

`nodes_analysis/08_Angular_Momentum_Transport_Node.md` records the derivation, the table,
the three verdicts, and the boundary: the balance law is checked but not integrated, there
is no gas model, no metric, and the radial part of the viscous divergence is not computed.

Mutation evidence (each reverted; both files byte-identical afterwards):

| Injected defect | Result |
| :--- | :--- |
| connection term `2 T / r` dropped from the balance | exit 1 |
| orbital shear off by a factor of two | exit 1 (2 tests) |
| orbital shear sign flipped | exit 1 (2 tests) |
| viscous stress sign flipped | exit 1 (6 tests) |
| dissipation factor dropped | exit 1 |
| `power_budget` always returns `sustained` | exit 1 |
| decay rate understated by nine orders of magnitude | exit 1 |
| `chi` given an `r-phi` component | exit 1 (2 tests) |

## Two repairs to files another thread had left mid-edit (2026-10-30)

The suite did not collect at the start of this session; both defects were in working-tree
changes made outside this task, and both were repaired minimally rather than reverted:

| Problem | Repair |
| :--- | :--- |
| `tests/test_ngc6503_vacuum_nodes.py` was reduced to a body with no imports (`NameError: pytest` at collection), used `np.trapz` (removed in NumPy 2.0), and read the origin's 0 → positive step as a sign crossing | restored the import/`RTOL` header, used `np.trapezoid` when available, excluded steps at `r <= 0.5` from the crossing search |
| `review/QVEX_NGC6503_Review_Report_v2.md` gained YAML front matter and lost its H1, so the one-title guard failed | promoted the document's own heading back to `#` |
