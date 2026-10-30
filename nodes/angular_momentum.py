"""
nodes/angular_momentum.py
The ``phi-r`` shear channel: angular-momentum transport by the vacuum stress.

Why this module exists
----------------------
Node 07 computes ``f^phi``, the azimuthal part of the exchange current. That is a *body*
torque supplied from outside the disk. A different claim is often made: that the vacuum
stress tensor itself, through its ``Pi^{phi r}`` component, carries angular momentum
radially and can therefore hold a flat rotation curve together. The two are not the same
mechanism, and this module separates them.

Two facts decide the question, and both are checkable.

1. **``Pi = Delta chi`` has no ``r-phi`` component.** With ``chi^{mu nu} = s^mu s^nu -
   h^{mu nu}/3`` and a radial ``s``, every off-diagonal term of ``chi`` vanishes
   identically: ``chi^{r phi} = 0``. A longitudinal anisotropy cannot torque anything,
   whatever its magnitude. To obtain ``Pi^{r phi} != 0`` the stress must have a *shear*
   part, ``Pi^{mu nu} = -2 eta sigma^{mu nu}``, with the orbital shear of a differentially
   rotating flow. That is the object implemented here.

2. **A tangential force cannot supply centripetal support.** The rotation curve is set by
   the *radial* force balance; ``Pi^{r phi}`` enters the ``phi`` equation. Transporting
   angular momentum reshapes ``Omega(r)``, it does not hold the gas up. The radial term is
   the one Node 06/07 writes as ``F_vac^r``, driven by ``p_r`` and ``Pi^{rr}``.

The orbital shear
-----------------
In coordinate spherical coordinates, with ``u^mu = gamma (1, 0, 0, Omega)``,

    sigma_{r phi} = (1/2)[d_r u_phi + d_phi u_r] - Gamma^lambda_{r phi} u_lambda
                  = (1/2)[d_r (r^2 gamma Omega) - (2/r) r^2 gamma Omega]
                  = (1/2) r^2 d_r (gamma Omega),

so the orthonormal component (``e_phi^ = (1/r) d_phi``) is

    sigma_{r^ phi^} = (1/2) r d_r (gamma Omega)  ->  (1/2) r dOmega/dr  for gamma = 1.

For solid-body rotation this vanishes; for a flat curve ``Omega = v/r`` it is ``-v/(2r)``;
for Keplerian rotation it is ``-3 Omega / 4``. The suite checks that derivation against a
finite-difference evaluation of the coordinate-basis expression, so the formula is not
taken on trust.

The balance
-----------
The ``phi`` component of ``nabla_mu T^{mu nu} = 0`` with only ``r``-dependence is

    (1/r^2) d_r (r^2 T^r{}_phi) = f_phi,

so ``r^2 T^r{}_phi`` is the angular-momentum flux and its divergence is the torque density.
For a flat curve with viscous stress ``T^r{}_phi = eta v / r`` the torque density is
``eta v / r^2`` and the power it must be paid for is ``eta v^2 / r^2``. The viscous
dissipation is ``2 eta sigma^2 = eta v^2 / (2 r^2)``: half the torque's work rate. Neither
is free, which is the quantitative part of the argument.

Units
-----
length kpc, time Myr, `Omega` 1/Myr, `sigma` 1/Myr, stress Msun/(kpc Myr^2),
`eta` Msun/(kpc Myr), torque density Msun/(kpc Myr^3), power density Msun/(kpc Myr^3).
``c = 1`` is *not* used: velocities are kpc/Myr, and the rest-energy comparison below takes
`c` explicitly.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

__all__ = [
    "LIGHT_SPEED_KPC_PER_MYR",
    "orbital_shear",
    "viscous_stress",
    "angular_momentum_flux",
    "torque_density",
    "dissipation",
    "flat_curve_stress",
    "flat_curve_torque_density",
    "flat_curve_power",
    "torque_power",
    "power_budget",
    "decay_rate_needed",
]

#: ``c`` in kpc/Myr, from c = 299792.458 km/s and 1 kpc = 3.0856775814913673e16 km.
LIGHT_SPEED_KPC_PER_MYR = 0.306601

_ATOL = 1e-15


def _positive(name: str, value):
    """Require every entry to be finite and strictly positive; return it as an array."""
    array = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(array)) or np.any(array <= 0.0):
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
    return array if array.ndim else float(array)


def orbital_shear(radius: float, omega: float, domega_dr: float, gamma: float = 1.0,
                  dgamma_dr: float = 0.0) -> float:
    """``sigma_{r phi} = (1/2) r d_r(gamma Omega)`` in 1/Myr.

    Vanishes for solid-body rotation, ``-v/(2r)`` for a flat curve ``Omega = v/r`` and
    ``-3 Omega/4`` for Keplerian rotation. ``gamma`` and its radial derivative are kept so
    the same expression covers the relativistic flow; for galactic speeds both are 1 and 0.
    """
    radius = _positive("radius", radius)
    omega = np.asarray(omega, dtype=float)
    domega_dr = np.asarray(domega_dr, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    dgamma_dr = np.asarray(dgamma_dr, dtype=float)
    for name, value in (("omega", omega), ("domega_dr", domega_dr),
                        ("gamma", gamma), ("dgamma_dr", dgamma_dr)):
        if np.any(~np.isfinite(value)):
            raise ValueError(f"{name} must be finite")
    if np.any(gamma <= 0.0):
        raise ValueError("gamma must be positive")
    result = 0.5 * radius * (gamma * domega_dr + omega * dgamma_dr)
    return float(result) if np.ndim(result) == 0 else result


def viscous_stress(eta: float, shear: float) -> float:
    """``Pi^{r phi} = -2 eta sigma^{r phi}`` in Msun/(kpc Myr^2).

    The sign follows from the Israel-Stewart form of Node 06: the stress opposes the shear.
    For a falling rotation curve (``dOmega/dr < 0``) the shear is negative and the stress is
    positive, so angular momentum is carried *outward* — the direction that flattens the
    curve. ``eta < 0`` would reverse that and is rejected: it would pump energy in.
    """
    eta = float(eta)
    if not np.isfinite(eta) or eta < 0.0:
        raise ValueError("shear viscosity must be non-negative and finite")
    shear = float(shear)
    if not np.isfinite(shear):
        raise ValueError("shear must be finite")
    return -2.0 * eta * shear


def angular_momentum_flux(radius: float, stress: float) -> float:
    """``r^2 T^r{}_phi``: the angular-momentum flux through a cylinder of radius ``r``.

    Units Msun/Myr (per unit height and per radian of azimuth, which cancels in the
    divergence below).
    """
    radius = _positive("radius", radius)
    return radius**2 * float(stress)


def torque_density(radius: float, stress: float, dstress_dr: float) -> float:
    """``(1/r^2) d_r(r^2 T^r{}_phi)`` in Msun/(kpc Myr^3): the torque per unit volume.

    This is the ``phi`` component of ``nabla_mu T^{mu nu}`` for an axisymmetric, height
    integrated flow. The ``2 T / r`` piece is the spherical connection term
    (``Gamma^mu_{mu r} - Gamma^phi_{r phi} = 1/r``), not a geometric correction to be
    dropped.
    """
    radius = _positive("radius", radius)
    stress = float(stress)
    dstress_dr = float(dstress_dr)
    if not np.isfinite(stress) or not np.isfinite(dstress_dr):
        raise ValueError("stress and its derivative must be finite")
    return 2.0 * stress / radius + dstress_dr


def dissipation(eta: float, shear: float) -> float:
    """``2 eta sigma^2`` in Msun/(kpc Myr^3): entropy produced, energy removed from the flow.

    The same ``Pi = -2 eta sigma`` that transports angular momentum dissipates. This is why
    the transport cannot be free, and why a flat curve maintained by viscosity needs an
    external energy input rather than merely a redistribution.
    """
    eta = float(eta)
    if not np.isfinite(eta) or eta < 0.0:
        raise ValueError("shear viscosity must be non-negative and finite")
    return 2.0 * eta * float(shear) ** 2


# --- the flat-curve case ----------------------------------------------------
def flat_curve_stress(eta: float, v_flat: float, radius: float) -> float:
    """``T^r{}_phi = eta v / r`` for ``Omega = v/r``, in Msun/(kpc Myr^2).

    ``sigma = -v/(2r)`` and ``T = -2 eta sigma``.
    """
    radius = _positive("radius", radius)
    v_flat = float(v_flat)
    if not np.isfinite(v_flat):
        raise ValueError("v_flat must be finite")
    if v_flat < 0.0:
        raise ValueError("v_flat is a speed and must be non-negative")
    return viscous_stress(eta, -0.5 * v_flat / radius)


def flat_curve_torque_density(eta: float, v_flat: float, radius: float) -> float:
    """``(1/r^2) d_r(r^2 T) = eta v / r^2`` for a flat curve, in Msun/(kpc Myr^3)."""
    radius = _positive("radius", radius)
    return flat_curve_stress(eta, v_flat, radius) / radius


def flat_curve_power(eta: float, v_flat: float, radius: float) -> float:
    """``f_phi v_phi = eta v^2 / r^2``: power per unit volume the torque must be paid for."""
    radius = _positive("radius", radius)
    v_flat = float(v_flat)
    return flat_curve_torque_density(eta, v_flat, radius) * v_flat


def torque_power(radius: float, torque: float, omega: float) -> float:
    """``tau * v_phi`` with ``v_phi = r Omega``: the work rate of a torque density."""
    radius = _positive("radius", radius)
    return float(torque) * radius * float(omega)


def power_budget(vacuum_power: float, torque: float, radius: float, omega: float,
                 eta: float = 0.0, shear: float = 0.0) -> Dict[str, object]:
    """Compare what the flow needs with what the exchange supplies.

    ``vacuum_power`` is the energy the exchange delivers per unit volume per unit time —
    Node 07's ``Q`` multiplied by the relevant energy scale, in Msun/(kpc Myr^3). The torque
    costs ``torque_power`` and the viscosity dissipates ``dissipation``; ``sustained`` is
    true only when the supply covers both.

    A negative ``net`` means the flow spins down: the vacuum shear redistributes angular
    momentum but the disk pays for it out of its own energy. That is the honest reading of
    the ``Pi^{r phi}`` channel, and it is why "flat rotation without an energy source" is
    not what this mechanism says.
    """
    radius = _positive("radius", radius)
    asked = torque_power(radius, float(torque), omega)
    lost = dissipation(eta, shear)
    total = asked + lost
    net = float(vacuum_power) - total
    return {
        "torque": torque,
        "torque_power": asked,
        "dissipation": lost,
        "required_power": total,
        "net": net,
        "sustained": bool(net >= -_ATOL * max(1.0, abs(total))),
    }


def decay_rate_needed(power: float, rest_energy_density: float) -> Dict[str, float]:
    """How fast the vacuum would have to convert rest energy to supply ``power``.

    ``rest_energy_density`` is ``rho_v c^2`` in the same units as ``power`` times a time,
    i.e. Msun/(kpc Myr^2). Returns the fractional rate (1/Myr) and the e-folding time (Myr).
    A rate comparable to or larger than the orbital frequency ``v/r`` means the vacuum is a
    transient, not a persistent halo — which is the test this number exists for.
    """
    power = float(power)
    rest_energy_density = _positive("rest_energy_density", rest_energy_density)
    if not np.isfinite(power):
        raise ValueError("power must be finite")
    if power <= 0.0:
        return {"rate_per_myr": 0.0, "efolding_myr": float("inf")}
    rate = power / rest_energy_density
    return {"rate_per_myr": rate, "efolding_myr": 1.0 / rate}
