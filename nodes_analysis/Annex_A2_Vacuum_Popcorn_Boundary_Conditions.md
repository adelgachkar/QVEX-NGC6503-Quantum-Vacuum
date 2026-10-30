# Annex A2: Vacuum Popcorn Model Physical Constraints & Analytical Boundary Conditions

- **Document ID:** `QVEX-ANNEX-A-BC-01`
- **Project:** QVEX (Quantum Vacuum Extended Analysis)
- **Target Node:** `nodes_analysis/01_Vacuum_Popcorn_Density_Node.md`
- **Epistemic Status:** **Unverified / Analytical Hypothesis Boundary Condition**
- **Note:** *Assuming a claim as a given hypothesis does not imply or establish its actual physical properties.*

---

## 1. Purpose and Epistemic Scope

This annex defines the rigorous analytical constraints, boundary limits, and consistency requirements for the vacuum popcorn density formulation within the QVEX framework.

The equations and parameters stated herein are treated strictly as
**analytical boundary conditions** and
**hypothetical test constraints**.

 They do not constitute verified empirical truths until validated against independent astrophysical and cosmological observables (such as rotation curves, lensing profiles, and CMB constraints).

---

## 2. Foundational Causal and Field Constraints

To ensure physical consistency and prevent unconstrained divergence, any implementation of the vacuum popcorn model must satisfy four mandatory foundational constraints:

### 2.1. Local Conservation & Bianchi Identity

$$\nabla^\mu T_{\mu\nu}^{(\text{popcorn})} = 0$$

- **Requirement:** The effective stress-energy tensor derived from the vacuum popcorn density fluctuations must be covariantly conserved, or coupled explicitly to a conserved source current if treated as an effective fluid.

### 2.2. Causal Horizon & Microcausality

- **Requirement:** Field fluctuations and correlation lengths must respect local microcausality. Vacuum popcorn bubble nucleation or coherence scales ($\xi$) cannot exceed the cosmological or local gravitational horizon at radius $r$:
  $$\xi(r) \le R_{\text{horizon}}(r)$$

### 2.3. Scale Decoupling & Boundary Screening

- **Requirement:** At asymptotic boundaries ($r \to \infty$ or high-density galactic cores where standard baryonic/dark matter dominates), the popcorn contribution must smoothly vanish or scale below observable detection thresholds without inducing abrupt discontinuities:
  $$\lim_{r \to \infty} \rho_{\text{popcorn}}(r) = 0$$

### 2.4. Radial-Tangential Stress Anisotropy Constraints

- **Requirement:** In spherically symmetric or axisymmetric configurations, the radial pressure ($p_r$) and tangential pressure ($p_t$) anisotropy must remain bounded:
  $$\frac{|p_r - p_t|}{\rho c^2} \le \Gamma_{\max}$$
  where $\Gamma_{\max}$ is bounded by local stability criteria (preventing infinite shear or runaway collapse).

---

## 3. Mathematical Formulation & Parameter Bounds

This annex governs Node 01, so it uses Node 01's symbols. Earlier revisions of this
section wrote the same profile as $r_s$, $r_c$ and $\alpha$; that second naming is
withdrawn because two names for one parameter is how the two documents drifted apart.
The baseline vacuum popcorn density profile is defined as:

$$\rho_{\text{popcorn}}(r) = \rho_0 \left( 1 + \frac{r}{\lambda_{\text{coh}}} \right)^{-\gamma_{\text{vac}}} \exp\left( -\frac{r}{R_{\text{cutoff}}} \right)$$

with the correspondence to the withdrawn symbols: $r_s \equiv \lambda_{\text{coh}}$,
$r_c \equiv R_{\text{cutoff}}$, $\alpha \equiv \gamma_{\text{vac}}$.

### 3.1 Canonical unit of the density amplitude

Canonical unit: **`kg/m^3`**. Node 01 quotes the amplitude in `g/cm^3`; that is a
conversion, not a different quantity (`1 g/cm^3 = 1000 kg/m^3`). The same value in every
supported unit:

| Parameter | Symbol | kg/m^3 (canonical) | g/cm^3 | Msun/kpc^3 |
| :--- | :--- | ---: | ---: | ---: |
| Amplitude | rho_0 | 7.02e-27 | 7.02e-30 | 1.037e2 |

A bound written in one unit is a bound on the density itself, so the constraint below is
stated in the canonical unit and holds whatever unit a reader prefers.

### Parameter Constraints

1. **Base Density ($\rho_0$):**
   - Range: $0 \le \rho_0 \le 10^{-25} \, \text{kg/m}^3$ (i.e. $\le 10^{-28} \, \text{g/cm}^3$)
   - Constraint: Must be strictly non-negative.
   - Note: Node 01's $7.02 \times 10^{-27} \, \text{kg/m}^3$ sits inside this range. Read as
     a `g/cm^3` number it would be a thousand times smaller, which is the error the
     canonical unit exists to make visible.
2. **Coherence Length ($\lambda_{\text{coh}}$):**
   - Range: $\lambda_{\text{coh}} > 0$ (measured in kpc for galactic applications).
3. **Cutoff Radius ($R_{\text{cutoff}}$):**
   - Range: $R_{\text{cutoff}} > \lambda_{\text{coh}}$ or $R_{\text{cutoff}} \to \infty$
     (for power-law dominant regimes).
4. **Decay Exponent ($\gamma_{\text{vac}}$):**
   - Range: $0 \le \gamma_{\text{vac}} \le 4$.

---

## 4. Implementation & Verification Instructions for PowerShell

To integrate and verify this annex in your local working environment, execute the following commands in PowerShell:

```powershell
# Ensure the nodes_analysis directory exists
if (!(Test-Path -Path "nodes_analysis")) {
New-Item -ItemType Directory -Path "nodes_analysis"
}

# Run pytest to verify existing project tests
python -m pytest -v
```

---

## 5. Document provenance (2026-10-29)

This document was recovered from the tail of `Annex_A_Vacuum_Popcorn_Physical_Constraints.md`,
where it had been appended inside an A.5 bullet of a different annex (QVEX-ANNEX-A-POP-01).
The two are separated because they name different targets, different baselines and different
epistemic statuses.

Later revisions change the symbols of section 3 and add section 3.1 (2026-10-30): the
profile is stated with Node 01's `lambda_coh`, `R_cutoff` and `gamma_vac`, the withdrawn
names are mapped one-to-one in the text, and the density amplitude is given in the
canonical unit with its `g/cm^3` and `Msun/kpc^3` equivalents. No numeric value changed;
only its unit was made explicit.

Structural changes made during the split, and only these: the heading `Annex A:` became
`Annex A2:`, the duplicated identifier `QVEX-ANNEX-A-POP-01` became `QVEX-ANNEX-A-BC-01`, the
unterminated ```powershell fence was closed, and this note was added. No technical text was
added, removed or reworded. The recovered section 4 ends where the source file resumed the
other annex; that is an artefact of the recovery, not a statement that the annex ends there.
