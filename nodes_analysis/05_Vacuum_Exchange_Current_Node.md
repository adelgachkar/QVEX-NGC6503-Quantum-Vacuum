# Node 05: Causal Exchange Current

## Objective

Replace the clamped Node 02 friction guard with a conserved exchange current, so that the
drag law used downstream cannot change sign with radius. Node 02's expression is retained
verbatim elsewhere as the documented formula, but it is no longer the production law.

The formal requirement comes from Annex A, A.3.3.2: a model that claims energy or momentum
exchange must state the system boundary, the source and sink terms, and a conservation
balance, including any exchange current. This node supplies that object.

---

## Mathematical Formulation (Canonical)

The vacuum-matter exchange coefficient as a function of galactocentric radius is

$$\Gamma(r) = \Gamma_0 \ln\left(1 + \frac{r}{r_{\text{core}}}\right)\left(1 - \exp\left(-\frac{r}{\lambda_{\text{coh}}}\right)\right)$$

### Parameters

```text
Gamma_0      = 1.1
r_core       = 4.5 kpc
lambda_coh   = 12.0 kpc
```

`Gamma_0`, `r_core` and `lambda_coh` take the values of the Node 02 activation profile. The
component of Node 02 that is *not* carried over is the linear penalty `-0.12 max(0, r-4.5)`.

### Exchange Current and Force

The two fluids exchange four-momentum through

$$J^{\mu} = \Gamma(r)\left(u_{\text{v}}^{\mu} - u_{\text{m}}^{\mu}\right)$$

with `+J` entering the matter fluid and `-J` entering the vacuum fluid. In the
non-relativistic limit the force on the matter fluid is

$$F(r) = -\Gamma(r)\,\Delta u$$

and the rate at which the exchange moves energy out of the matter fluid is
`u_m · J = -Gamma · (du)^2`.

---

## Why the Sign Change Cannot Recur

For `r >= 0` both factors of `Gamma` are non-negative, and the coefficient is zero only at
`r = 0`. The sign of the force is therefore carried entirely by `Delta u`, the relative
velocity between the two fluids. Consequences, each of which is machine-checked:

1. `Gamma(r) >= 0` on the whole domain, so the exchange cannot act as a drive.
2. For a fixed relative velocity the force has one sign at every radius; it does not reverse
   near 12.11 kpc or anywhere else.
3. The dissipated power is `-Gamma (du)^2 <= 0`, so energy only leaves the matter sector.
4. The matter source and the vacuum sink sum to zero identically, not approximately.

The failure mode of Node 02 was not a calibration problem: the penalty made the expression
itself negative, and a clamp then hid the fact. Here positivity is a property of the
definition, so no guard is required and none is present.

---

## Domain and Validation

`r` must be finite and non-negative. Negative radii are rejected with `ValueError` rather
than clipped, following Annex A, A.3.1.2 (do not silently reinterpret invalid input).

## What This Node Does Not Claim

- `Gamma(r)` is a calibrated profile, not a derived one: no field equation produces it, and
  its shape is inherited from Node 02 rather than predicted.
- Dropping the penalty changes the outer-disk magnitude. Any previous fit that relied on
  the penalised expression must be re-run; the numbers in `ngc6503_metrics.json` are
  residuals against a fixed baryonic reference and are unaffected because no mapping from
  friction to circular velocity exists yet.- The exchange current is a conservation structure. It is not evidence for vacuum energy
transfer, and it does not by itself produce a rotation-curve prediction.
- **The two scalars above are a truncation of one four-vector, and the dropped part is not
  always small.** `F = -Gamma du` and `u_m.J = -Gamma (du)^2` are the `O(dv)` and `O(dv^2)`
  pieces of `J^mu = Q u^mu + f^mu`; the azimuthal component `f^phi`, which carries torque,
  is absent. At a front with `v_r = 0.05` and `v_phi = 0.15` (vacuum at rest, `Gamma = 1`)
  the current is `Q = +0.0127`, `f^r = -0.0513`, `f^phi = -0.1539`: the part this node
  reports is a third of the momentum transfer and none of the torque. See
  `07_Exchange_Four_Vector_Conservation_Node.md`, which keeps all four components and
  reduces to the formulas above as `dv -> 0`.

### Extension: Quantum-Hydrodynamic Compressive Barrier and Effective Potential

## Section X: Vacuum Popcorn Density Law and Exchange Current Framework

### X.1 Definitions and Field Parameters

The baseline vacuum popcorn model characterizes the spatial distribution of the vacuum energy density and the associated exchange current parameters within the QVEX-NGC6503 framework. The fundamental parameters are defined as:

- $\rho_0$: Base vacuum density scale ($\rho_0 \approx 7.02 \times 10^{-30} \text{ g/cm}^3$).
- $r_{\text{core}}$: Core characteristic radius ($4.5 \text{ kpc}$).
- $\lambda_{\text{coh}}$: Coherence length ($\lambda_{\text{coh}} = 12.0 \text{ kpc}$).
- $R_{\text{cutoff}}$: Outer boundary cutoff radius ($R_{\text{cutoff}} = 50.0 \text{ kpc}$).
- $\gamma_{\text{vac}}$: Density falloff exponent ($\gamma_{\text{vac}} = 1.0$).

### X.2 Vacuum Popcorn Density Distribution

The radial density profile of the vacuum popcorn medium is modeled via a modified power-law decay modulated by an exponential cutoff:
$$\rho_{\text{popcorn}}(r) = \rho_0 \left(1 + \frac{r}{\lambda_{\text{coh}}}\right)^{-\gamma_{\text{vac}}} \exp\left(-\frac{r}{R_{\text{cutoff}}}\right)$$
*Constraint & Validity Domain:* Defined strictly for $r \ge 0$. Negative radii trigger a `ValueError` per Annex A boundary protocols.

### X.3 Vacuum Exchange Current Law (Node 05 Formulation)

To resolve the unphysical negative friction regime encountered in earlier iterations (Node 02), the exchange current damping function $\Gamma(r)$ is formulated as:
$$\Gamma(r) = \Gamma_0 \ln\left(1 + \frac{r}{r_{\text{core}}}\right)\left(1 - \exp\left(-\frac{r}{\lambda_{\text{coh}}}\right)\right)$$
where $\Gamma_0 = 1.1$ represents the exchange amplitude. The former piecewise clamp $\max(0, r - r_{\text{core}})$ is officially deprecated and removed.

### X.4 Stress-Energy Tensor and Pressure Constraints (Status: Auxiliary / Non-Operational)

- **Conservation Law:** The effective vacuum exchange satisfies the covariant conservation constraint $\nabla^\mu T_{\mu\nu}^{(\text{popcorn})} = 0$.
- **Stress Anisotropy Bound:** Radial and tangential pressures are bounded by the stress-anisotropy limit:
  $$\frac{|p_r - p_t|}{\rho c^2} \le \Gamma_{\text{max}}$$
- **Unresolved / Pending Derivations:** A rigorous derivation of a dynamic bulk modulus ($\mathcal{B}$), compressive barrier potential ($V_{\text{comp}}$), and direct rotational velocity mapping ($v_c(r)$) remains unmapped in the current v0.1.0 codebase and requires explicit tensor-level coupling beyond calibrated phenomenological friction.
def compressive_acceleration(self, r: np.ndarray, v_sound: float = 120.0) -> np.ndarray:
    """Computes radial outward compressive acceleration a_comp(r) in (km/s)^2 / kpc."""
    r_arr = np.asarray(r, dtype=float)
    denom = self.r_core**2 + r_arr**2
    return (self.alpha *v_sound**2)* r_arr) / denom

def circular_velocity_qvex(
    self,
    r: np.ndarray,
    v_bar: np.ndarray,
    v_sound: float = 120.0,
    gamma_scale: float = 1.0,
) -> np.ndarray:
    """Calculates total rotation curve v_c(r) = sqrt(v_bar^2 + v_vac^2 + v_comp^2 + v_gamma^2)."""
    r_arr = np.asarray(r, dtype=float)
    v_vac_sq = (self.G *self.effective_enclosed_mass(r_arr)) / np.maximum(r_arr, 1e-6)
    v_comp_sq = r_arr* self.compressive_acceleration(r_arr, v_sound=v_sound)
    gamma_r = self.exchange_coefficient(r_arr)
    v_gamma_sq = gamma_scale *gamma_r* r_arr
    v_total_sq = np.maximum(0.0, v_bar**2 + v_vac_sq + v_comp_sq + v_gamma_sq)
    return np.sqrt(v_total_sq)
