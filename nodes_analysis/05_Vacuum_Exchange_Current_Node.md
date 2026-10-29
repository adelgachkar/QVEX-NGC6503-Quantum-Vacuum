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
  friction to circular velocity exists yet.
- The exchange current is a conservation structure. It is not evidence for vacuum energy
  transfer, and it does not by itself produce a rotation-curve prediction.
