# Node 08: Angular-Momentum Transport by the `phi-r` Shear

- **Document ID:** `QVEX-NODE-08-TORQUE-01`
- **Project:** QVEX (Quantum Vacuum Extended Analysis)
- **Implements:** `nodes/angular_momentum.py`
- **Bound by:** `tests/test_angular_momentum.py`
- **Depends on:** Node 06 (`Pi^{mu nu}`, transport), Node 07 (`f^phi`, conservation)
- **Epistemic Status:** **Unverified / Effective Description** — the algebra and the
  balance law are exact for the stated stress; the rotation-curve numbers are order-of-
  magnitude estimates, not a fit.
- **Note:** *Assuming a claim as a given hypothesis does not imply or establish its actual physical properties.*

---

## 1. The claim under test

The proposition usually stated is this: if the leading and trailing fronts produce a
non-zero shear stress `Pi^{phi r}`, angular momentum is transferred directly from the
active vacuum to the rotating gas; if `(1/sqrt(-g)) d_r(sqrt(-g) Pi^{phi r})` balances the
friction and the pressure gradient, `Omega_m(r)` can hold a flat rotation curve without
Keplerian fall-off and without dark matter; and the four-dimensional formulation constrains
not only the radial pressure deficit but the rotation itself.

Three separable assertions hide in that sentence, and they do not stand or fall together:

| # | Assertion | Verdict |
| :--- | :--- | :--- |
| A | A `phi-r` vacuum stress transports angular momentum | **True, and now implemented** — but not with the tensor currently used |
| B | That transport can hold a flat rotation curve together | **False as stated** — a tangential force supplies no centripetal support |
| C | Therefore no dark matter (or no energy source) is needed | **False** — the torque must be paid for, and the price is measurable |

## 2. Assertion A: the channel exists, but not in `Pi = Delta chi`

With `chi^{mu nu} = s^mu s^nu - h^{mu nu}/3` and a radial `s`, every off-diagonal component
of `chi` is **zero**, identically and for every `Delta`:

$$s = (0,1,0,0) \;\Longrightarrow\; \chi^{r\phi} = s^r s^\phi - \tfrac13 h^{r\phi} = 0 - 0 = 0$$

So the `Pi^{mu nu}` of Nodes 06 and 07 — the one whose `Pi : sigma` enters the energy
identity — cannot torque anything. A longitudinal anisotropy is a *normal* stress; a torque
needs a *shear* stress. The suite asserts the zero exactly (`tensor[1, 3] == 0.0`), not to
within a tolerance, because no parameter choice changes it.

To obtain `Pi^{r phi} != 0` the stress must carry the Israel-Stewart shear part,

$$\Pi^{\mu\nu} = -2\eta_{\text{vac}}\,\sigma^{\mu\nu}, \qquad \sigma^{\mu\nu} = \tfrac12 h^{\mu\alpha}h^{\nu\beta}\left(\nabla_\alpha u_\beta + \nabla_\beta u_\alpha\right) - \tfrac13 \Theta h^{\mu\nu},$$

which is non-zero only for **differential** rotation. That is the object this node adds.

## 3. The orbital shear, derived rather than quoted

In coordinate spherical coordinates with `u^mu = gamma (1, 0, 0, Omega)` and
`u_phi = r^2 gamma Omega`:

$$\sigma_{r\phi} = \tfrac12\left[\partial_r u_\phi + \partial_\phi u_r\right] - \Gamma^\lambda{}_{r\phi}u_\lambda = \tfrac12\left[\partial_r(r^2\gamma\Omega) - \tfrac{2}{r}r^2\gamma\Omega\right] = \tfrac12 r^2 \partial_r(\gamma\Omega)$$

and the orthonormal component is `sigma = (1/2) r d_r(gamma Omega)`, i.e. for galactic
speeds

$$\sigma_{r\hat\phi} = \tfrac12 r \frac{d\Omega}{dr}.$$

The `2/r` term is the connection `Gamma^lambda_{r phi} u_lambda`, not a correction to be
dropped; without it the result is off by a factor of three. The suite checks this against a
finite-difference evaluation of the same coordinate expression, so the formula is derived in
the test rather than trusted from the module.

Consequences, each tested:

| Rotation | `sigma_{r phi}` | `T^{r phi} = -2 eta sigma` |
| :--- | :--- | :--- |
| Solid body (`Omega` const) | `0` | `0` — the channel needs differential rotation |
| Flat curve (`Omega = v/r`) | `-v/(2r)` | `+eta v/r` — outward transport |
| Keplerian (`Omega ∝ r^{-3/2}`) | `-3 Omega/4` | `+3 eta Omega/2` |

The sign is the physically correct one: for a falling rotation curve the stress is positive,
so angular momentum flows **outward**, which is the direction that flattens `Omega(r)`.
Reversing `dOmega/dr` reverses the transport.

## 4. Assertion B: a tangential force cannot hold the gas up

The `phi` component of `nabla_mu T^{mu nu} = 0` for an axisymmetric flow is

$$\frac{1}{r^2}\frac{d}{dr}\left(r^2 T^r{}_\phi\right) = f_\phi,$$

so `T^r{}_phi` enters the **azimuthal** equation and nothing else. The radial equation — the
one that sets `v_c` — contains `T^{rr}`, `p_r` and `Pi^{rr}`. Transporting angular momentum
reshapes the profile `Omega(r)`; it does not supply centripetal acceleration `v^2/r`.

This is why "flatten the curve with `Pi^{phi r}`" cannot work as an explanation of the mass
deficit: the deficit is radial, and a torque is tangential. What the torque *can* do is
govern how the profile evolves — spin-up, spin-down, or the viscous spreading of a ring.

One caveat, stated so that it is not mistaken for a rescue: the divergence of the shear
stress is *not* purely azimuthal, because the same `Pi` has a `phi phi` component and the
frame connection `omega^{r}{}_{phi phi} = -1/r` converts it into a radial term. That term is
`O(eta Omega / r)` — second order in the flow rather than a pressure gradient — and this node
does **not** implement or bound it. It is recorded as a limitation, not used as an argument
in either direction.

## 5. Assertion C: the price of the torque

For a flat curve the stress is `T^r{}_phi = eta v/r`, so

$$\text{torque density} = \frac{\eta v}{r^2}, \qquad P_{\text{torque}} = \tau\,\Omega = \frac{\eta v^2}{r^2}, \qquad \mathcal{D} = 2\eta\sigma^2 = \frac{\eta v^2}{2r^2}.$$

`P_torque = 2 D` — the classical disk relation, checked rather than quoted. Both must be
supplied from outside the flow; neither is free.

Numbers for NGC 6503 (`v_flat = 115 km/s = 0.1175 kpc/Myr`), the canonical vacuum amplitude
`rho_v = 103.7 Msun/kpc^3` (Node 01) and the largest viscosity the causality bound of Node 06
allows, `eta = tau (rho_v c^2) = 195 Msun/(kpc Myr)` at `tau = 20 Myr` (`c = 0.3066 kpc/Myr`,
so `rho_v c^2 = 9.75 Msun/(kpc Myr^2)`):

| `r` (kpc) | `sigma` (1/Myr) | `T^{r phi}` | torque | `P_torque` | `D` | total | rate (1/Myr) | e-folding (Myr) | orbit (Myr) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | -1.175e-02 | 4.582 | 0.9163 | 0.1077 | 0.0538 | 0.1615 | 1.657e-02 | 60 | 267 |
| 10 | -5.875e-03 | 2.291 | 0.2291 | 0.0269 | 0.0135 | 0.0404 | 4.142e-03 | 241 | 535 |
| 20 | -2.938e-03 | 1.145 | 0.0573 | 0.0067 | 0.0034 | 0.0101 | 1.036e-03 | 966 | 1070 |

(units: stress and `rho c^2` Msun/(kpc Myr^2); torque Msun/(kpc Myr^3); power
Msun/(kpc Myr^3); rate 1/Myr.)

Two readings follow, and both matter:

1. **The vacuum would have to decay on a timescale shorter than one rotation.** At `r = 10`
   kpc the e-folding time of 241 Myr is less than the orbital period of 535 Myr. A halo that
   converts a third of its rest energy per orbit is not a static halo. Even at the causality
   ceiling the mechanism is a **transient**, so it cannot be the steady explanation of a
   rotation curve that has looked the same for gigayears.
2. **The torque is small where the numbers are honest.** At `r = 10` kpc the shear stress is
   `2.29 Msun/(kpc Myr^2)`, which is `6.6e-3`, `1.7e-3` and `1.7e-4` of the gas ram pressure
   `rho_gas v^2` for `rho_gas = 2.5e4`, `1e5` and `1e6 Msun/kpc^3`. To reach 10% of the ram
   pressure at `1e5 Msun/kpc^3` the viscosity would have to be `1.2e4`, sixty times the
   causality ceiling — i.e. excluded.

The scaling is the part worth keeping: `Q_required ∝ eta`, and so is the torque. A smaller
`eta` is not cheaper in a way that helps — it makes the mechanism weaker by exactly the same
factor. There is no regime in which the transport is both dynamically relevant and
energetically cheap.

## 6. What is implemented, and what is not

| Element | Function | Property tested |
| :--- | :--- | :--- |
| Orbital shear | `orbital_shear` | matches the coordinate-basis covariant derivative; zero for solid body; `-v/2r` flat, `-3 Omega/4` Keplerian |
| Viscous stress | `viscous_stress` | `-2 eta sigma`; sign reverses with `dOmega/dr`; `eta < 0` raises |
| Angular-momentum flux | `angular_momentum_flux` | `r^2 T^r{}_phi`, positive (outward) for a falling curve |
| Balance law | `torque_density` | equals the finite-difference divergence of the flux |
| Dissipation | `dissipation` | `2 eta sigma^2 >= 0`; `P_torque = 2 D` for a flat curve |
| Flat-curve case | `flat_curve_stress`, `flat_curve_torque_density`, `flat_curve_power` | closed forms consistent with the balance law |
| Energy ledger | `power_budget` | reports the deficit instead of assuming support; exact break-even |
| Cost in decay rate | `decay_rate_needed` | rate scales with `eta`; e-folding versus orbital period |

**Not implemented — the boundary of this node:**

- **No evolution.** The balance law is checked, not integrated: nothing here follows
  `Omega(r, t)` forward in time, so "can it *maintain* a flat curve" is answered by an
  energy budget, not by a simulation.
- **No gas physics.** The disk has no surface density profile, no star-formation or cooling
  term, and `rho_gas` appears only in the comparison above as a stated number.
- **The radial part of the viscous divergence is not computed.** See section 4: the `phi phi`
  component of the same stress feeds a radial term through the connection, and that term is
  outside what this module evaluates.
- **No field-equation coupling.** `Pi^{mu nu}` is prescribed through `eta` and `Omega(r)`;
  `G_{mu nu} = 8 pi G T_{mu nu}/c^4` is not solved, so no metric follows from the torque.
- **`eta_vac` is not derived.** Node 06's causality bound is used as an upper limit; where in
  that range the vacuum actually sits is a free parameter, and the table is computed at the
  limit that flatters the mechanism most.
- **No observational comparison.** `v_flat = 115 km/s` and `r = 5, 10, 20 kpc` are read off
  the SPARC curve for orientation, not fitted.

## 7. Units

length kpc, time Myr, `Omega` and `sigma` 1/Myr, stress and `rho c^2` Msun/(kpc Myr^2),
`eta` Msun/(kpc Myr), torque and power densities Msun/(kpc Myr^3). `c = 0.306601 kpc/Myr`
explicitly (`LIGHT_SPEED_KPC_PER_MYR`); unlike Nodes 06 and 07 this node does not set
`c = 1`, because its central number is a ratio of an energy rate to a rest-energy density.

## 8. Relation to the earlier nodes

Node 06 defined `Pi` but only its longitudinal part, and recorded that the angular-momentum
sector was missing. Node 07 kept `f^phi`, which is a *body* torque from the exchange, not a
stress. This node supplies the stress channel and separates the two: `f^phi` pushes the gas
from outside, `Pi^{r phi}` redistributes what the gas already has. Neither provides the
radial support that a flat rotation curve actually requires.
