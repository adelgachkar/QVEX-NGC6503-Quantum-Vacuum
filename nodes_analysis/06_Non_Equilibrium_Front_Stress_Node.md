# Node 06: Non-Equilibrium Stress of Leading and Trailing Fronts

- **Document ID:** `QVEX-NODE-06-FRONT-01`
- **Project:** QVEX (Quantum Vacuum Extended Analysis)
- **Implements:** `nodes/non_equilibrium_stress.py`
- **Bound by:** `tests/test_non_equilibrium_front_stress.py`
- **Epistemic Status:** **Unverified / Effective Description** — the algebra below is
  exact for the stated assumptions; no parameter is fitted to an observation.
- **Note:** *Assuming a claim as a given hypothesis does not imply or establish its actual physical properties.*

---

## 1. The gap this node closes

Nodes 01, 02 and 05 describe an **equilibrium** vacuum: a static density profile
`rho(r)`, a static exchange coefficient `Gamma(r)`, and a drag force
`F = -Gamma * dv` that answers instantly at every radius. Such a law is a function of
position only. It therefore cannot distinguish a front that is *arriving* from one that
has *passed*: both sit at the same `r` and get the same `Gamma`.

A front is exactly the place where that assumption fails. The vacuum state is changing
there, so the stress it carries lags the shear that produces it. The lag is not a
correction of order one percent — as the table in section 4 shows, it can remove 78% of
the stress the equilibrium law would predict, and it is the difference between a
compression and a rarefaction.

## 2. Geometry: what makes `Pi^{mu nu}` well defined

Work in an orthonormal tetrad `(e_t, e_r, e_theta, e_phi)` with signature `(-,+,+,+)`, so
that no metric factors enter the contractions:

$$u^\mu u_\mu = -1, \qquad s^\mu s_\mu = +1, \qquad u^\mu s_\mu = 0$$

`u^mu` is the four-velocity of the matter, `s^mu` the unit radial direction of the local
rest frame. The spatial projector and the traceless longitudinal projector are

$$h^{\mu\nu} = g^{\mu\nu} + u^\mu u^\nu = \mathrm{diag}(0,1,1,1)$$

$$\chi^{\mu\nu} = s^\mu s^\nu - \frac{1}{3} h^{\mu\nu}$$

so that the dissipative stress of an anisotropic vacuum is

$$\Pi^{\mu\nu} = \Delta \, \chi^{\mu\nu}, \qquad \Delta \equiv p_r - p_t$$

All three defining properties then hold *identically*, not approximately:

$$u_\mu \Pi^{\mu\nu} = 0, \qquad \Pi^\mu{}_\mu = \Delta\left(1 - \tfrac{3}{3}\right) = 0, \qquad \Pi^{\mu\nu} = \Pi^{\nu\mu}$$

In components `Pi_rr = 2 Delta / 3` and `Pi_theta-theta = Pi_phi-phi = -Delta / 3`. The
inverse relation used by the code is `Delta = (3/2) Pi_rr`.

This is a stronger statement than the one usually made in prose. "Traceless" and
"orthogonal to `u`" are not conditions imposed on a computed tensor and then checked; they
are consequences of the basis, which is why `stress_tensor()` cannot return a tensor that
violates them.

## 3. Transport: the stress has memory

Israel-Stewart relaxation for the shear part, in the fluid frame:

$$\tau_\Pi \, h^\mu{}_\alpha h^\nu{}_\beta \, u^\rho \nabla_\rho \Pi^{\alpha\beta} + \Pi^{\mu\nu} = -2\eta_{\text{vac}} \, \sigma^{\mu\nu}$$

with `tau_Pi` the relaxation time of the vacuum fluctuations, `eta_vac` the effective
shear viscosity and `sigma^{mu nu}` the shear of the flow. For a fluid element crossing a
front of width `w` at speed `v0`, the crossing time is `t_cross = w / v0` and, starting
from `Pi = 0` with a constant shear, the exact solution is

$$\Pi(t) = -2\eta_{\text{vac}}\sigma\left(1 - e^{-t/\tau_\Pi}\right)$$

$$\Pi(t_{\text{cross}}) = -2\eta_{\text{vac}}\sigma\left(1 - e^{-\mathrm{De}}\right), \qquad \mathrm{De} \equiv \frac{\tau_\Pi}{t_{\text{cross}}}$$

`De` is the Deborah number. `De << 1`: the stress equilibrates before the front passes and
the equilibrium law is right. `De >= 1`: the stress has not arrived; the front is
genuinely non-equilibrium and the stress is deposited *behind* it.

Causality is a constraint, not a remark: the shear signal speed
`v_c = sqrt(eta / (tau_Pi (rho + p)))` must satisfy `v_c <= 1`. Taking `tau_Pi -> 0`
recovers the parabolic law `Pi = -2 eta sigma`, which propagates instantly and is therefore
excluded. The implementation refuses such parameters instead of silently integrating them.

**Sign of the dissipation.** With the dissipative stress entering as
`T^{mu nu} = ... - 2 eta sigma^{mu nu}`, the contraction `Pi : sigma` is *negative* for a
dissipative front. The second law requires `-(1/2) Pi : sigma = eta sigma : sigma >= 0`.
Reporting `Pi : sigma` itself as entropy production would claim the front runs backwards in
time; `entropy_production()` returns the correctly signed quantity and the suite tests it.

## 4. Two fronts from one velocity field

An outward-moving shell has two fronts. With `r_in < r_out`,

$$v(r) = \frac{v_0}{2}\left[\tanh\frac{r - r_{\text{in}}}{w} - \tanh\frac{r - r_{\text{out}}}{w}\right]$$

which is `0` outside, `v0` inside, and smooth across both edges. In the thin-front,
non-relativistic limit `sigma_rr = (2/3) dv/dr`, so at each front's own centre
`|sigma| = v0 / (3 w)` with

$$\sigma_{\text{leading}} = -\frac{v_0}{3w} < 0, \qquad \sigma_{\text{trailing}} = +\frac{v_0}{3w} > 0$$

The leading (outer) front decelerates the flow into the ambient medium — a compression. The
trailing (inner) front accelerates it away from rest — a rarefaction. The two fronts
therefore carry **opposite** dissipative stress. This is the statement the Node 01/02/05
chain cannot make, because `Gamma(r)` has one sign at every `r`.

Numbers for `w = 0.5 kpc`, `v0 = 0.1 kpc/Myr` (`~98 km/s`), `eta = 0.3`, `rho + p = 10`
(in the node's units) — leading front only, `t_cross = 5 Myr`:

| `tau_Pi` (Myr) | `De` | `sigma_rr` | `Pi_rr` equilibrium | `Pi_rr` at crossing | reached | `-(1/2)Pi:sigma` | `v_c` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.2 | 0.04 | -0.0667 | 0.0400 | 0.0400 | 100.0% | 1.33e-3 | 0.387 |
| 2.0 | 0.40 | -0.0667 | 0.0400 | 0.0367 | 91.8% | 1.22e-3 | 0.123 |
| 20.0 | 4.00 | -0.0667 | 0.0400 | 0.0088 | 22.1% | 2.95e-4 | 0.039 |

Reproduce with `front_summary(10.0, 20.0, 0.5, 0.1, eta=0.3, tau=tau, inertia=10.0)`.

## 5. The interaction engine

Total conservation is an identity (Bianchi), so the matter source is fixed by the vacuum
sector rather than chosen:

$$\nabla_\mu T_{\text{tot}}^{\mu\nu} = 0 \;\Longrightarrow\; \nabla_\mu T_m^{\mu\nu} = -\nabla_\mu T_{\text{vac}}^{\mu\nu} = Q^\nu$$

$$Q^\nu = Q\,u^\nu + f^\nu, \qquad u_\nu f^\nu = 0$$

$$Q = -u_\nu \nabla_\mu T_m^{\mu\nu}, \qquad f^\alpha = h^\alpha{}_\nu Q^\nu$$

and the projected matter equation is the generalized Euler equation

$$(\rho_m c^2 + p_m) a^\alpha + h^{\alpha\mu}\nabla_\mu p_m = f^\alpha, \qquad a^\alpha = u^\mu \nabla_\mu u^\alpha$$

whose source is the divergence of the anisotropic vacuum stress,

$$F_{\text{vac}}^\alpha \equiv -h^\alpha{}_\nu \nabla_\mu\left[(p_r - p_t)s^\mu s^\nu + \Pi^{\mu\nu}\right]$$

Two consequences worth stating explicitly:

1. If the vacuum is isotropic and uniform, `grad_mu T_vac^{mu nu} ∝ grad^nu Lambda = 0` and
   the bracket vanishes: the vacuum applies no local force, whatever its energy density.
2. With `Pi^{mu nu} != 0` the bracket is non-zero even where `p_r - p_t = 0`, because
   `Pi` itself is a gradient source. That is the term the equilibrium nodes approximate
   away, and the term that can torque the flow.

## 6. What is implemented, and what is not

Implemented in `nodes/non_equilibrium_stress.py`, and bound by the suite:

| Element | Function | Property tested |
| :--- | :--- | :--- |
| Tetrad, normalization, orthogonality | `geometry`, `four_velocity`, `spatial_direction` | `u.u = -1`, `s.s = +1`, `u.s = 0`; a tilted `s` is rejected |
| Projectors, `chi`, `Pi = Delta chi` | `projector`, `anisotropy_tensor`, `stress_tensor` | `u.Pi = 0`, `tr Pi = 0`, symmetry, component split |
| Equilibrium limit | `equilibrium_stress` | `Pi = -2 eta sigma`, antisymmetric in `sigma` |
| Israel-Stewart solution | `relaxed_stress` | matches the analytic exponential at `t = tau` and `t -> inf` |
| Front kinematics | `shell_velocity`, `radial_shear`, `front_shear` | opposite signs, sampled peak equals the analytic peak |
| Deborah number, crossing time | `crossing_time`, `deborah_number` | `De << 1` -> equilibrium, `De >> 1` -> non-equilibrium; monotone |
| Causality | `causal_speed`, `validate_transport` | `v_c > 1`, `tau <= 0`, `eta < 0` all raise |
| Second law | `entropy_production` | `-(1/2)Pi:sigma >= 0` at both fronts |
| Interaction vector | `exchange_decomposition` | `u.Q = -Q`, `f` orthogonal to `u`, matter + vacuum = 0 |

**Not implemented — the honest boundary of this node:**

- **No equations of state.** `p_r = -rho c^2 + Delta_r` and `p_t = -rho c^2 + Delta_t` are
  parameters of the model, not derived from a vacuum microphysics.
- **No metric and no field-equation solve.** `sigma^{mu nu}` is taken from a prescribed
  radial velocity field, not from `G_{mu nu} = 8 pi G (T_m + T_vac)/c^4`. Solving that
  system requires a metric ansatz (Kerr-like or with `g_{t phi}` for rotation) and a
  boundary-value solve that this bundle does not contain.
- **No angular momentum sector in *this module*.** `Sigma` (shear) and `Omega` (vorticity)
  are dropped and the tetrad here is radial only, so the torque that `f^phi` would exert is
  outside it. Node 07 (`nodes/exchange_current.py`) keeps the azimuthal component of the
  current and computes `r f^phi`; Node 08 (`nodes/angular_momentum.py`) adds the shear part
  `-2 eta sigma` whose `phi r` component is the one that actually transports angular
  momentum — note that the `Pi = Delta chi` defined here has `Pi^{r phi} = 0` identically.
  Neither later node supplies a metric or a field-equation solve.
- **No observational binding.** Nothing here is compared with a rotation curve. The
  parameters `eta`, `tau_Pi`, `w`, `v0` are free, and the numbers in section 4 are a
  demonstration of the mechanism, not a fit.
- **`T_m^{mu nu}` is not specified.** `rho_m` and `p_m` appear only as the inertia in the
  causality bound.

## 7. Units

length kpc, time Myr, velocity kpc/Myr, stress and pressure Msun/(kpc Myr^2), shear
viscosity `eta` Msun/(kpc Myr), relaxation time `tau_Pi` Myr, shear `sigma` 1/Myr.
`c = 1` in the causality bound. The vacuum density itself stays in the canonical unit of
`01_Vacuum_Popcorn_Density_Node.md` (`kg/m^3`, converted explicitly); this node never
re-derives it.

## 8. Relation to the earlier nodes

Node 05's exchange current remains the production drag law for a *steady* flow. Node 06
adds the part Node 05 lacks: a front is not steady, and the difference between the
equilibrium stress and the delivered stress is `1 - exp(-De)`. For `De << 1` the two nodes
agree, which is the regime in which Node 05's drag law is meaningful; for `De >= 1` any
prediction built on the instantaneous law overstates the stress, at one front and
understates it at the other.
