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

The baseline vacuum popcorn density profile is defined as:

$$\rho_{\text{popcorn}}(r) = \rho_0 \left( 1 + \frac{r}{r_s} \right)^{-\alpha} \exp\left( -\frac{r}{r_c} \right)$$

### Parameter Constraints

1. **Base Density ($\rho_0$):**
   - Range: $0 \le \rho_0 \le 10^{-25} \, \text{kg/m}^3$
   - Constraint: Must be strictly non-negative.
2. **Scale Radius ($r_s$):**
   - Range: $r_s > 0$ (typically measured in kpc for galactic applications).
3. **Core/Cutoff Radius ($r_c$):**
   - Range: $r_c > r_s$ or $r_c \to \infty$ (for power-law dominant regimes).
4. **Decay Exponent ($\alpha$):**
   - Range: $0 \le \alpha \le 4$.

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

Structural changes made during the split, and only these: the heading `Annex A:` became
`Annex A2:`, the duplicated identifier `QVEX-ANNEX-A-POP-01` became `QVEX-ANNEX-A-BC-01`, the
unterminated ```powershell fence was closed, and this note was added. No technical text was
added, removed or reworded. The recovered section 4 ends where the source file resumed the
other annex; that is an artefact of the recovery, not a statement that the annex ends there.
