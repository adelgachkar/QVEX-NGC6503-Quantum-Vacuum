# Node 07: The Exchange Four-Vector and Total Conservation

- **Document ID:** `QVEX-NODE-07-EXCHANGE-01`
- **Project:** QVEX (Quantum Vacuum Extended Analysis)
- **Implements:** `nodes/exchange_current.py`
- **Bound by:** `tests/test_exchange_current.py`
- **Depends on:** Node 05 (`J^mu = Gamma (u_v - u_m)`), Node 06 (`Pi^{mu nu}`, transport)
- **Epistemic Status:** **Unverified / Effective Description** — the identities below are
  exact for the stated tensor; no parameter is fitted to an observation.
- **Note:** *Assuming a claim as a given hypothesis does not imply or establish its actual physical properties.*

---

## 1. The gap: a four-vector with one component

Node 05 defines the exchange current correctly and then keeps two numbers from it,

$$J^\mu = \Gamma\left(u_{\text{v}}^\mu - u_{\text{m}}^\mu\right) \;\longrightarrow\; F = -\Gamma\,\Delta u, \qquad u_{\text{m}}\!\cdot\!J = -\Gamma(\Delta u)^2,$$

a drag force and an energy rate. Both are **radial and scalar**. The azimuthal component
`f^phi` — the one that carries torque — is not approximated; it is dropped before it is
ever computed. Node 06 then notes in its own limitations section that "the tetrad is radial
only, so the torque that `f^phi` would exert is outside the model."

This node closes that gap without adding a parameter. The current keeps all four
components, the split into energy and force is enforced rather than assumed, and the
torque becomes a computed number.

## 2. Frame conditions

In an orthonormal tetrad `(e_t, e_r, e_theta, e_phi)` with signature `(-,+,+,+)`:

$$u^\mu u_\mu = -1, \qquad s^\mu s_\mu = +1, \qquad u^\mu s_\mu = 0, \qquad h^{\mu\nu} = g^{\mu\nu} + u^\mu u^\nu, \qquad \chi^{\mu\nu} = s^\mu s^\nu - \tfrac{1}{3}h^{\mu\nu}$$

The radial axis of the *rest frame of `u`* is not `e_r` when the flow moves radially. It is

$$s^\mu = \frac{e_r^\mu + (u\cdot e_r)\,u^\mu}{\sqrt{1 + (u\cdot e_r)^2}},$$

which satisfies `s.s = +1` and `u.s = 0` identically, including the time component that a
radially moving flow forces on it. Handing `e_r` to such a flow leaves a rest-frame energy
flux in an *equilibrium* configuration; the suite rejects it instead of accepting it
(`validate_frame` raises on `u.s != 0`).

**An index-position trap, and the bug it caused here.** The first version of this module
projected a contravariant vector with the all-up tensor,

$$f^\alpha \stackrel{?}{=} h^{\alpha\nu}J_\nu \quad\text{applied as}\quad h^{\mu\nu}J^\nu,$$

which is not the same object: `(g + u\!\otimes\!u) u != 0`. The projector that acts on a
vector has one index up and one down,

$$h^\alpha{}_\nu = \delta^\alpha{}_\nu + u^\alpha u_\nu, \qquad f^\alpha = J^\alpha + u^\alpha (u\cdot J),$$

and `(h^{mu nu} u_nu) = 0` is its covector twin. The two coincide only in the rest frame
of `u`, which is exactly why a rest-frame-only implementation hides the error. Both are
now separate functions (`spatial_projector` and `flow_projector`), and the suite checks
both.

## 3. The current

$$J^\mu = Q\,u^\mu + f^\mu, \qquad u_\mu f^\mu = 0, \qquad Q = -u_\mu J^\mu, \qquad f^\alpha = h^\alpha{}_\nu J^\nu$$

`Q > 0` means the vacuum sector loses energy to the matter sector. `ExchangeCurrent`:

- refuses a force with a component along `u` (`f.u != 0` raises) — orthogonality is a
  definition here, not a check applied afterwards;
- splits *any* four-vector exactly (`ExchangeCurrent.from_vector`), so nothing is lost on the
  way in;
- reports `torque_density(r) = r f^phi`, which is zero for a radial force by construction.

**Reduction of Node 05.** Expanding `J = Gamma(u_v - u_m)` for a small relative velocity:

| `dv` | `Q` | `Q / (Gamma dv^2 / 2)` | `f^r / (-Gamma dv)` |
| ---: | ---: | ---: | ---: |
| 1e-2 | 5.000e-05 | 1.000075006249 | 1.000100010001 |
| 1e-3 | 5.000e-07 | 1.000000750239 | 1.000001000001 |
| 1e-4 | 5.000e-09 | 0.999999976669 | 1.000000010000 |

So Node 05's two scalars are the `O(dv)` and `O(dv^2)` pieces of one current, and the
remainder is not small in general: for a rotating flow the azimuthal force can dominate.

## 4. Total conservation fixes the source

`nabla_mu T_tot^{mu nu} = 0` follows from the Bianchi identity, so the matter source is not
a choice:

$$\nabla_\mu T_{\text{m}}^{\mu\nu} = -\nabla_\mu T_{\text{vac}}^{\mu\nu} = J^\nu$$

For the anisotropic, dissipative vacuum stress

$$T_{\text{vac}}^{\mu\nu} = (\rho_{\text{v}} + p_t)u^\mu u^\nu + p_t g^{\mu\nu} + \Delta\, s^\mu s^\nu + \Pi^{\mu\nu}, \qquad \Pi^{\mu\nu} = \Pi_0\,\chi^{\mu\nu},$$

contracting the divergence with `u` is an identity. Term by term:

| Piece of `u_nu d_mu T^{mu nu}` | Contribution |
| :--- | :--- |
| `u_nu d_mu[(rho + p_t) u^mu u^nu]` | `-[rho_dot + p_t_dot + (rho + p_t) Theta]` |
| `u_nu d_mu[p_t g^{mu nu}]` | `+p_t_dot` |
| `u_nu d_mu[Delta s^mu s^nu]` | `-Delta (s^mu s^nu d_mu u_nu)` |
| `u_nu d_mu[Pi^{mu nu}]` | `-Pi : sigma` |

and the `p_t_dot` cancels, leaving

$$\boxed{\;Q = \dot\rho_{\text{v}} + (\rho_{\text{v}} + p_t)\,\Theta + \Delta\,\left(s^\mu s^\nu \nabla_\mu u_\nu\right) + \Pi : \sigma\;}$$

with `rho_dot = u^mu d_mu rho` the **convective** derivative, `Theta = nabla.u` the
expansion, and `Pi : sigma` the dissipative contraction.

Three things are worth naming, because each is a place where a simpler model silently
differs:

1. **The pressure in the `p dV` term is `p_t`, not the mean.** An anisotropic vacuum does
   work differently along its two axes; using `(p_r + 2p_t)/3` changes the answer even when
   `Delta = 0`.
2. **The `Delta` term is the stretching of the anisotropy axis along itself.** A purely
   radial model sets `s^mu s^nu d_mu u_nu` to a fixed function of `r` and cannot see that
   the axis tilts as the flow rotates.
3. **`Pi : sigma` survives where `Delta = 0`.** The dissipative stress is its own source, so
   a non-equilibrium front exchanges energy even in a perfectly isotropic vacuum.

## 5. The projected matter equation

$$(\rho_{\text{m}} + p_{\text{m}})a^\alpha + h^{\alpha\mu}\nabla_\mu p_{\text{m}} = f^\alpha, \qquad a^\alpha = u^\mu\nabla_\mu u^\alpha$$

$$a^\alpha = -\frac{1}{\rho_{\text{m}} + p_{\text{m}}}\left[h^{\alpha\mu}\nabla_\mu p_{\text{m}} - F_{\text{vac}}^\alpha\right], \qquad F_{\text{vac}}^\alpha = -h^\alpha{}_\nu \nabla_\mu\left[(p_r - p_t)s^\mu s^\nu + \Pi^{\mu\nu}\right]$$

`euler_acceleration()` returns that `a^alpha` from the force and the projected pressure
gradient, and refuses a force that is not orthogonal to the flow. Two limits are tested:

- `F_vac = 0`: `a = -h grad p / (rho + p)`, orthogonal to `u` — ordinary pressure support.
- `grad p = 0`: `a = f / (rho + p)` — the vacuum force is the whole acceleration.

If the vacuum is isotropic and uniform, `grad_mu T_vac^{mu nu} ∝ grad^nu Lambda = 0` and the
bracket vanishes: no local force, whatever the energy density.

## 6. The torque sector

`f^phi` is what turns momentum transfer into angular momentum transfer. At a front where the
matter moves with `v_r = 0.05` and `v_phi = 0.15` while the vacuum is at rest, with
`Gamma = 1` (`w = 0.5 kpc`, `v0 = 0.1 kpc/Myr`, `r_in = 10`, `r_out = 20`):

| `r` (kpc) | `v_r` | `Q` | `f^r` | `f^phi` | torque `r f^phi` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10.0 | 0.0500 | +0.01274 | -0.05128 | -0.15385 | -1.5385 |
| 20.0 | 0.0500 | +0.01274 | -0.05128 | -0.15385 | -3.0769 |

The same force at two radii gives two torques, as it must. More to the point: the azimuthal
component is **three times** the radial one here, and Node 05's scalar law reports only the
radial part. Reversing the drift (vacuum ahead of the matter) flips the sign of the torque
and the two cases differ by exactly the boost factor `gamma` of the rotating flow, because
they are the same current seen from two rest frames.

What this does **not** yet give is a rotation curve. `a^phi != 0` is a torque density, not
a `v_c(r)`: closing that requires a metric with `g_{t phi}`, a boundary-value solve, and a
stellar/gas tracer model.

Nor is `f^phi` the same thing as a stress: it pushes the gas from outside, while a
`Pi^{phi r}` component would redistribute the gas's own angular momentum. That second
channel is Node 08, and it is *not* available from the `Pi = Delta chi` used here — see the
section on `chi^{r phi}` there.

## 7. Numerical evidence that the identity is the one implemented

The identity of section 4 is checked against an independent finite-difference divergence of
`vacuum_stress` on a `(t, r)` grid, with a front that moves outward and a flow that carries
an azimuthal component:

| grid `n x n` | relative residual | signal `|Q|` |
| ---: | ---: | ---: |
| 41 x 41 | 1.50e-02 | 0.1569 |
| 81 x 81 | 4.82e-03 | 0.1625 |
| 161 x 161 | 2.75e-03 | 0.1642 |

The residual falls with refinement while the signal stays fixed, so what is left is
discretization, not a missing term. The test asserts `< 2e-2` at `81 x 81`; the measured
value is `4.8e-3`.

Mutations, each of which turns the suite red (and each reverted):

| Injected defect | Result |
| :--- | :--- |
| `Pi` dropped from the stress tensor | 1 failed |
| `p_t` replaced by the isotropic mean | 1 failed |
| projection built from the all-up tensor | 6 failed |
| torque sign reversed | 2 failed |
| `d_t` dropped from the divergence | 1 failed |
| `s = e_r` handed to a moving flow | 11 failed |

## 8. What is implemented, and what is not

| Element | Function | Property tested |
| :--- | :--- | :--- |
| Frame conditions | `validate_frame`, `boosted_velocity`, `radial_direction` | `u.u = -1`, `s.s = +1`, `u.s = 0`; superluminal and tilted frames raise |
| Projectors | `spatial_projector`, `flow_projector`, `project_onto_flow_orthogonal` | `h^mu{}_nu u^nu = 0`, `h^{mu nu} u_nu = 0`, `tr h = 3` |
| Current and split | `ExchangeCurrent` and its classmethod `ExchangeCurrent.from_vector` | `f.u = 0`, `u.J = -Q`, exact round trip |
| Reduction of Node 05 | `exchange_from_relative_velocity` | `Q = Gamma dv^2/2` and `f^r = -Gamma dv` at small `dv` |
| Torque | `torque_density`, `rotating_front_current` | sign flips with the drift; radial force gives zero |
| Conservation ledger | `conservation_ledger` | legs sum to zero; equal legs do not balance |
| Euler equation | `euler_acceleration` | orthogonal to `u`; both limits |
| Energy identity | `vacuum_energy_exchange`, `flat_divergence` | matches the numerical divergence to discretization order |
| Rest-frame energy | `vacuum_stress_energy_in_frame` | `u.T = -rho u` exactly, for any `Delta`, `Pi` |

**Not implemented — the boundary of this node:**

- **No equations of state.** `p_r`, `p_t`, `rho_v` and `Pi_0` are parameters, not derived
  from vacuum microphysics. The `Delta_r`, `Delta_t` decomposition is a bookkeeping choice.
- **No metric and no field-equation solve.** The divergence is taken in a flat orthonormal
  frame, with `d_theta = d_phi = 0` by axisymmetry. `G_{mu nu} = 8 pi G T_{mu nu}/c^4` is not
  solved, so nothing here is a solution of Einstein's equations.
- **No closure of the matter sector.** `rho_m` and `p_m` appear only as the inertia in the
  Euler equation; `T_m` is not specified, and the system `(4) + (5)` is not solved as a
  boundary-value problem.
- **No observational binding.** `eta`, `tau_Pi`, `Gamma`, `w`, `v0` are free. The numbers in
  section 6 are a demonstration that the sector exists and has the right sign, not a fit.
- **No microcausality check on the exchange itself.** `J` is not required to be timelike or
  null; only the Israel-Stewart signal speed of Node 06 is bounded.

## 9. Units

length kpc, time Myr, velocity kpc/Myr, `T^{mu nu}` Msun/(kpc Myr^2), `J^mu`
Msun/(kpc^2 Myr^2) (a stress per length, i.e. a force density), torque density `r f^phi`
Msun/(kpc Myr^2), inertia `rho + p` Msun/(kpc Myr^2). `c = 1`. The vacuum density itself
stays in the canonical unit of `01_Vacuum_Popcorn_Density_Node.md`.

## 10. Relation to the earlier nodes

Node 05 remains the equilibrium drag law, and this node agrees with it in the limit it was
derived in (`dv -> 0`, radial flow). Node 06 supplies the non-equilibrium `Pi` that enters
the energy identity through `Pi : sigma`. What is new here is the *vector* structure: the
split `Q u + f` with `f` orthogonal to `u`, the energy identity with `p_t` and both
anisotropy terms, the Euler equation that turns the vacuum divergence into an acceleration,
and the azimuthal component that makes the torque computable at all.
