# Annex A — Vacuum-Popcorn Physical Constraints

**Document ID:** QVEX-ANNEX-A-POP-01  
**Project:** QVEX — Quantum Vacuum Exchange Engine  
**Baseline:** QVEX v0.1.0  
**Document class:** Model-governance and physical-consistency annex  
**Epistemic status:** Class C phenomenological working model; not experimentally validated  
**Language:** English  
**Purpose:** Normative constraints for implementation, calibration, reporting, and interpretation of the vacuum-popcorn model.

---

## A.1. Status and scope

This annex defines the minimum physical, dimensional, numerical, and evidentiary conditions for using the vacuum-popcorn formulation in QVEX analyses. “Normative” means required for a result to be described as compliant with this project model; it does **not** mean that the underlying vacuum-nucleation hypothesis is established physics.

The baseline computational node is `qvex.nodes.popcorn_vacuum_density.PopcornVacuumDensityNode`. It accepts baryonic acceleration `g_bar` and calculates reference acceleration scales and a phenomenological compensation profile. It does not, by itself, simulate bubble nucleation, solve a field equation, establish a stress-energy tensor, or demonstrate a causal vacuum interaction.

> **Bundle provenance (2026-10-29).** The baseline module named above,
> `qvex.nodes.popcorn_vacuum_density.PopcornVacuumDensityNode`, is **not shipped** in the
> `QVEX-NGC6503-Quantum-Vacuum` bundle; it lives in the parent engine QVEX v0.1.0
> (`QVEX_v0.1.0/qvex/nodes/popcorn_vacuum_density.py`). The bundle's executable node is
> `nodes.active_tensor_foam.ActiveTensorFoamSolver`, whose profile is governed by
> `nodes_analysis/01_Vacuum_Popcorn_Density_Node.md` and Annex A2. Therefore the statements
> in A.2 and A.4.3 that `lambda_coh_kpc` and `vacuum_energy_density` are inactive describe the
> parent engine only: in this bundle `lambda_coh` and `rho_0` are live parameters of
> `vacuum_density()`. Neither set of parameters acts on the other package.


These requirements apply to model calculations and their reporting. They do not replace general relativity, quantum field theory, established observational constraints, or laboratory safety and metrology requirements. The separate axial-jet/dipole hypothesis and PID parameter mapping are optional extensions, not part of the baseline popcorn node, unless a study explicitly implements and validates those couplings.

## A.2. Definitions and baseline relations

Use the following notation consistently:

- `g_bar`: baryonic acceleration, in m s⁻²; it must be finite and non-negative.
- `H0`: Hubble parameter supplied in km s⁻¹ Mpc⁻¹.
- `a0`: reference acceleration scale computed from `H0`.
- `g_pop`: model popcorn acceleration contribution, in m s⁻².
- `f_nuc`: dimensionless vacuum-nucleation/compensation factor used by the code. It is a model diagnostic, not a measured nucleation rate or probability.
- `lambda_coh_kpc`: configured coherence-length parameter in kpc. The v0.1.0 node stores this parameter but does not use it in the calculation shown below.

The baseline implementation uses:

\[
H_{0,\mathrm{SI}} = H_0\frac{1000\ \mathrm{m\,s^{-1}}}{1\ \mathrm{Mpc}},\qquad
 a_0(H_0)=\frac{c H_{0,\mathrm{SI}}}{2\pi},
\]

\[
g_{\mathrm{pop}}=\sqrt{g_{\mathrm{bar}}a_0},\qquad
f_{\mathrm{nuc}}=\sqrt{\frac{a_{0,\mathrm{local}}}{\max(g_{\mathrm{bar}},10^{-15}\ \mathrm{m\,s^{-2}})}}.
\]

The code defaults are `H0_early = 67.4 km s⁻¹ Mpc⁻¹` and `H0_late = 73.04 km s⁻¹ Mpc⁻¹`. They are reference inputs, not immutable constants. The configured `lambda_coh_kpc` default is 12.0 kpc. The dataclass declares `vacuum_energy_density = 5.35e-10`, but the v0.1.0 computation does not use that field and does not specify its unit in the node implementation. It must not be treated as a measured or operative energy density without an explicit unit definition, provenance, and implemented coupling.

## A.3. Mandatory physical and dimensional constraints

### A.3.1. Units and input domain

1. Convert `H0` from km s⁻¹ Mpc⁻¹ to s⁻¹ before calculating `a0`; use the exact conversion `1 Mpc = 3.085677581491367 × 10²² m` and `c = 299,792,458 m s⁻¹` or documented, consistent precision.
2. `g_bar` must be supplied in m s⁻². Values expressed in another unit must be converted before entering the node. Reject NaN, infinity, and negative accelerations; do not silently reinterpret them.
3. `a0` has units m s⁻². The radicand `g_bar × a0` has units m² s⁻⁴, so `g_pop` has units m s⁻². The ratio in `f_nuc` is dimensionless.
4. Any additional term combined with `g_bar` or `g_pop` must have acceleration units. An energy density, pressure, rate, or dimensionless factor cannot be added to an acceleration without a stated, dimensionally valid coupling.
5. Preserve the unit and provenance of every fitted parameter in configuration, output metadata, plots, and tables. Do not infer units merely from a variable name.

### A.3.2. Positivity and limiting behavior

1. For the baseline model, `a0 > 0` when `H0 > 0`; `g_pop ≥ 0` for every accepted `g_bar ≥ 0`.
2. `g_pop` must approach zero as `g_bar` approaches zero, and the baseline formula gives `g_pop = 0` at exactly zero.
3. At fixed positive `a0`, `g_pop` is monotonic in `g_bar` and obeys `g_pop² = g_bar × a0`. The factor `f_nuc` is not itself an acceleration or a directly observed quantity.
4. The numerical floor `10⁻¹⁵ m s⁻²` in the denominator of `f_nuc` is a regularization convention. It prevents division by zero; it does not establish a physical lower bound on acceleration, a minimum nucleation rate, or a vacuum-density floor. Results sensitive to the floor must be identified and tested against alternative floors or an explicit zero-input policy.
5. Do not extrapolate the phenomenological square-root relation beyond the calibrated or explicitly studied domain without labeling the result as extrapolation.

### A.3.3. Conservation, gravity, and interpretation

1. A positive `g_pop` is an effective model contribution. It is not, without an independently specified field equation and source term, evidence for a new force, mass component, vacuum-energy transfer, or gravitationally active bubble population.
2. If a model claims energy or momentum exchange, it must state the system boundary, source and sink terms, and a conservation balance. For a covariant relativistic model, specify the stress-energy tensor and the applicable conservation equation, including any exchange current; do not imply conservation from a plotted acceleration curve.
3. If `g_total = g_bar + g_pop` is reported, state the sign convention, geometry, and whether the terms are assumed to be collinear. Use a vector or potential-based calculation when that assumption is not justified.
4. No force, metric perturbation, thrust, or causal vacuum-exchange claim follows from the baseline popcorn node alone.

## A.4. Parameter and implementation controls

1. Record the software version, source revision, configuration, input-data identifier, units, and numerical precision for each run. Keep defaults distinct from fitted values and observational measurements.
2. A parameter may be called “calibrated” only when the calibration dataset, objective function, priors or bounds, uncertainty treatment, and validation split are documented. Parameters adjusted on the evaluation data are not independent validation.
3. `vacuum_energy_density` and `lambda_coh_kpc` are presently declared configuration values but are not used by the baseline calculation. Report them as inactive unless the exact source version demonstrates otherwise. Do not imply that varying an inactive parameter changes the prediction.
4. The node compares early- and local-reference `H0` values through `a0` and `g_pop`; this is a parameterized comparison, not a derivation or resolution of the Hubble tension.
5. The optional NGC 6503 axial dipole/jet potential, nonlinear foam-friction function, and PID gains documented elsewhere in the research notes are not automatically coupled to this node. If combined, document each equation, dimensional coefficient, boundary condition, and data path; avoid treating controller gains as physical constants.
6. Changes to equations, defaults, floors, or output names require a versioned change record and corresponding regression tests.

## A.5. Numerical and data-quality requirements

- Use finite floating-point inputs and verify output finiteness. Preserve array shape and input ordering.
- Test scalar and vector inputs, the zero limit, values around the regularization floor, the high-acceleration regime, and the low-acceleration regime.
- Check monotonicity and the identity `g_pop**2 ≈ g_bar * a0` within a tolerance justified by numerical precision.
- Propagate uncertainties in `H0`, baryonic acceleration, distances, inclinations, mass-to-light ratios, and relevant calibration parameters. Report covariance where material; do not present a best-fit curve without uncertainty.
- For galaxy rotation-curve work, use documented observational data and quality selections. Separate training/calibration galaxies or points from held-out validation data; account for distance and inclination systematics and compare with appropriate baseline models.
- Do not count software unit tests, synthetic examples, or a successful optimizer run as observational evidence. Tests establish implementation behavior only.

## A.6. Validation, falsification, and reporting

Before interpreting a result physically, the study must:

1. State the exact equations and distinguish implemented terms from proposed mechanisms.
2. Provide a dimensional audit and limiting-case tests.
3. Compare predictions with independent data and at least one relevant competing or null model using a declared metric and uncertainty model.
4. Pre-register, or clearly mark as post hoc, the observable and rejection threshold. A discrepancy criterion must include the full uncertainty budget and treatment of correlated errors.
5. Report failed fits, parameter degeneracies, sensitivity to data selection, and sensitivity to numerical regularization.
6. Use language such as “phenomenological model contribution,” “hypothesis,” or “fit” unless independent evidence supports a stronger statement.

The hypothesis is not supported by a good fit alone. It must yield reproducible, discriminating predictions and survive independent tests. A result inconsistent with a pre-specified prediction at the declared confidence level counts against that tested model version; it must not be rescued by undocumented parameter changes.

## A.7. Compliance checklist

A result is compliant with this annex only when all applicable items below are satisfied:

- [ ] English-only technical record with equation and software-version identifiers.
- [ ] Inputs and outputs have explicit units; invalid or non-finite inputs are rejected.
- [ ] Dimensional, positivity, limiting-case, and regularization checks pass.
- [ ] Inactive parameters and phenomenological assumptions are disclosed.
- [ ] Fit, calibration, and independent validation data are distinguished.
- [ ] Uncertainties, systematic effects, null model, and falsification criterion are reported.
- [ ] Claims do not exceed what the code and evidence demonstrate.

## A.8. Version note

This annex describes the QVEX v0.1.0 popcorn-vacuum node as inspected for this document. The stated numerical defaults reflect that implementation and may change in later releases. This document is a project constraint and reporting standard, not a claim that vacuum-popcorn nucleation has been experimentally observed or that QVEX is a validated physical engine.


---

## A.9. Document provenance (2026-10-29)

This file previously carried a second, unrelated document appended inside the A.5 bullet list.
That document is now `Annex_A2_Vacuum_Popcorn_Boundary_Conditions.md` with identifier
`QVEX-ANNEX-A-BC-01`. The split is structural: no technical text was added, removed or
reworded. Changes made, and only these: the appended document was moved out (restoring the
A.5 list), the bundle-provenance note in A.1 was added, and this record was added.
