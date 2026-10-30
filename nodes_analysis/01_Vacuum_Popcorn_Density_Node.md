# Node 01: Vacuum Popcorn Density and Zero-Point Fluctuations

## Objective

Define the spatial distribution of the quantum vacuum popcorn energy density in the outer boundaries of galactic halos and establish its preliminary parameterization.

---

## Mathematical Formulation (Canonical)

The popcorn vacuum energy density profile as a function of galactocentric radius $r$ is defined as:

$$\rho_{\text{popcorn}}(r) = \rho_0 \left(1 + \frac{r}{\lambda_{\text{coh}}}\right)^{-\gamma_{\text{vac}}} \exp\left(-\frac{r}{R_{\text{cutoff}}}\right)$$

### Parameters

```text
rho_0        = 7.02 × 10^-30 g/cm^3
lambda_coh   = 12.0 kpc
gamma_vac    = 1.0
R_cutoff     = 50.0 kpc
```

### Canonical Unit of the Density Amplitude

The profile above is dimensionless except for `rho_0`, which carries the whole mass
dimension. The amplitude is therefore stated once, in one canonical unit, and every other
representation is a conversion of it — never a second definition.

Canonical unit: **`kg/m^3`**. The `g/cm^3` value in the parameter block is the same
quantity expressed in the historical unit, and the two are related exactly by
`1 g/cm^3 = 1000 kg/m^3`.

| Parameter | Symbol | kg/m^3 (canonical) | g/cm^3 | Msun/kpc^3 |
| :--- | :--- | ---: | ---: | ---: |
| Amplitude | rho_0 | 7.02e-27 | 7.02e-30 | 1.037e2 |

The two length scales, `lambda_coh` and `R_cutoff`, are in `kpc`; the profile itself is
unit-free in `r/lambda_coh` and `r/R_cutoff`, so no further conversion enters the
equation. A consumer that needs `Msun/kpc^3` (for example to combine the density with
`G = 4.30091e-6 kpc (km/s)^2/Msun`) must convert first:

```text
rho_0 [Msun/kpc^3] = rho_0 [g/cm^3] * 1000 / 6.7682e-29
                   = 7.02e-30 * 1.4775e31 = 1.037e2
```

`Annex_A2_Vacuum_Popcorn_Boundary_Conditions.md` states the same amplitude in `kg/m^3`;
`tests/test_density_units.py` compares the two documents after conversion, so a number
placed in the wrong unit column turns the suite red instead of passing as a plausible
looking value.
