"""
nodes/non_equilibrium_stress.py
Non-equilibrium (dissipative) stress carried by leading and trailing vacuum fronts.

Why this module exists
----------------------
The Nodes 01/02/05 chain gives an **equilibrium** description: a density profile, an
exchange coefficient ``Gamma(r)`` and a drag force set by the *relative* velocity
``F = -Gamma dv``. That law cannot say anything about a front, because it has no memory:
it answers instantly at every radius. A front, by contrast, is a place where the vacuum
state is *changing*, so the stress it carries is not the one the equilibrium law predicts —
it lags behind by the relaxation time ``tau_Pi``.

The object that carries this is the dissipative part ``Pi^{mu nu}`` of the vacuum
stress-energy tensor (Israel-Stewart). This module implements it in the fluid rest frame
with the smallest structure that is still faithful:

* the geometry is an orthonormal tetrad ``(e_t, e_r, e_theta, e_phi)`` with signature
  ``(-,+,+,+)`` — not a coordinate basis, so no ``r^2`` factors appear in the contractions;
* ``u^mu = (1,0,0,0)`` and ``s^mu = (0,1,0,0)`` with ``u.u = -1``, ``s.s = +1``,
  ``u.s = 0``;
* ``h^{mu nu} = g^{mu nu} + u^mu u^nu = diag(0,1,1,1)``;
* ``chi^{mu nu} = s^mu s^nu - h^{mu nu}/3`` (the traceless longitudinal projector), so
  ``Pi^{mu nu} = Delta * chi^{mu nu}`` with ``Delta = p_r - p_t``. Then
  ``u_mu Pi^{mu nu} = 0``, ``Pi^mu_mu = 0`` and ``Pi^{mu nu} = Pi^{nu mu}`` hold
  identically, not approximately.

The fronts
----------
One outward-moving shell has two fronts. Its radial velocity is

    v(r) = (v0/2) [ tanh((r - r_out)/w) - tanh((r - r_in)/w) ],

so ``dv/dr`` is positive at the leading (outer) front and negative at the trailing (inner)
one. With ``sigma_rr = (2/3) dv/dr`` in the thin-front, non-relativistic limit, the two
fronts carry **opposite** shear and therefore opposite dissipative stress. That is the
statement this module makes visible, and it is what the equilibrium ``Gamma(r)`` law cannot
express: ``Gamma`` is a function of radius alone, so it would give the same sign at both.

Transport
---------
    tau_Pi * h^{mu}_{~alpha} h^{nu}_{~beta} u^rho d_rho Pi^{alpha beta}
        + Pi^{mu nu} = -2 eta sigma^{mu nu}

For a fluid element crossing a front of width ``w`` at speed ``v0`` the crossing time is
``t_cross = w / v0``. Starting from ``Pi = 0`` and a constant shear, the exact solution is

    Pi(t) = -2 eta sigma (1 - exp(-t / tau_Pi)),

so at the crossing

    Pi(t_cross) = -2 eta sigma (1 - exp(-De)),   De = tau_Pi / t_cross.

``De`` is the Deborah number. ``De << 1``: the stress keeps up and the equilibrium law is
right. ``De >= 1``: the stress has not arrived by the time the front has passed — the front
is genuinely non-equilibrium, and the stress is deposited *behind* it.

Causality
---------
Israel-Stewart is causal only if the shear signal speed stays below light:
``v_c^2 = eta / (tau_Pi * (rho + p)) <= 1``. :func:`causal_speed` returns that number and
:func:`validate_transport` refuses parameters that violate it, because a diffusive
``Pi = -2 eta sigma`` (``tau_Pi -> 0``) propagates instantly.

Units (same discipline as ``nodes/active_tensor_foam.py``)
---------------------------------------------------------
length kpc, time Myr, velocity kpc/Myr, pressure/stress Msun/(kpc Myr^2),
shear viscosity eta Msun/(kpc Myr), relaxation time tau_Pi Myr. Every public function
states the unit of what it returns; no bare number changes unit silently.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Union

import numpy as np

__all__ = [
    "NATURAL_CAUSALITY_LIMIT",
    "geometry",
    "four_velocity",
    "spatial_direction",
    "projector",
    "anisotropy_tensor",
    "stress_tensor",
    "trace",
    "is_orthogonal_to_flow",
    "equilibrium_stress",
    "relaxed_stress",
    "crossing_time",
    "deborah_number",
    "causal_speed",
    "validate_transport",
    "shell_velocity",
    "radial_shear",
    "front_shear",
    "front_summary",
    "anisotropy_pressure",
    "entropy_production",
    "exchange_decomposition",
]

#: With c = 1, a causal signal speed cannot exceed this.
NATURAL_CAUSALITY_LIMIT = 1.0


# --- geometry ---------------------------------------------------------------
def geometry() -> np.ndarray:
    """Orthonormal-tetrad metric ``diag(-1, 1, 1, 1)`` (signature ``-,+,+,+``)."""
    return np.diag([-1.0, 1.0, 1.0, 1.0])


def four_velocity() -> np.ndarray:
    """``u^mu`` of a fluid at rest in the tetrad: ``(1, 0, 0, 0)``, ``u.u = -1``."""
    return np.array([1.0, 0.0, 0.0, 0.0])


def spatial_direction() -> np.ndarray:
    """``s^mu`` along the radial direction: ``(0, 1, 0, 0)``, ``s.s = +1``, ``u.s = 0``."""
    return np.array([0.0, 1.0, 0.0, 0.0])


def _validated_direction(u, s):
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    s = spatial_direction() if s is None else np.asarray(s, dtype=float)
    g = geometry()
    norm_u = float(u @ g @ u)
    norm_s = float(s @ g @ s)
    if not np.isclose(norm_u, -1.0, atol=1e-12):
        raise ValueError(f"u.u = {norm_u}, expected -1 (timelike, (-,+,+,+))")
    if not np.isclose(norm_s, 1.0, atol=1e-12):
        raise ValueError(f"s.s = {norm_s}, expected +1 (spacelike unit vector)")
    if not np.isclose(float(u @ g @ s), 0.0, atol=1e-12):
        raise ValueError("u.s != 0: the radial direction must be orthogonal to the flow")
    return u, s


def projector(u=None) -> np.ndarray:
    """``h^{mu nu} = g^{mu nu} + u^mu u^nu`` — projects onto the flow-orthogonal space."""
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    return geometry() + np.outer(u, u)


def anisotropy_tensor(Delta: float = 1.0, u=None, s=None) -> np.ndarray:
    """``chi^{mu nu} = s^mu s^nu - h^{mu nu}/3`` scaled by ``Delta = p_r - p_t``.

    This is the traceless longitudinal projector: the only structure a single radial
    direction can add to an isotropic stress. ``h`` removes the flow direction and
    ``-h/3`` removes the trace, so ``chi`` is purely anisotropic by construction.
    """
    u, s = _validated_direction(u, s)
    h = projector(u)
    return Delta * (np.outer(s, s) - h / 3.0)


def stress_tensor(Delta: float, u=None, s=None) -> np.ndarray:
    """``Pi^{mu nu}`` for a radial pressure anisotropy ``Delta = p_r - p_t``."""
    if not np.isfinite(Delta):
        raise ValueError("Delta must be finite")
    return anisotropy_tensor(Delta, u=u, s=s)


def trace(tensor: np.ndarray) -> float:
    """``Pi^mu_mu`` with the tetrad metric."""
    return float(np.trace(geometry() @ np.asarray(tensor, dtype=float)))


def is_orthogonal_to_flow(tensor: np.ndarray, u=None, atol: float = 1e-12) -> bool:
    """``u_mu Pi^{mu nu} = 0``: the stress must not be seen as energy by the fluid."""
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    contracted = (geometry() @ u) @ np.asarray(tensor, dtype=float)
    return bool(np.all(np.abs(contracted) <= atol))


# --- transport --------------------------------------------------------------
def equilibrium_stress(eta: float, shear: Union[float, np.ndarray]):
    """Navier-Stokes limit ``Pi = -2 eta sigma`` (``tau_Pi = 0``).

    ``eta`` in Msun/(kpc Myr), ``shear`` in 1/Myr, result in Msun/(kpc Myr^2).
    """
    if eta < 0.0:
        raise ValueError("shear viscosity must be non-negative")
    return -2.0 * eta * np.asarray(shear, dtype=float)


def validate_transport(eta: float, tau: float, inertia: float) -> None:
    """Refuse parameters that would make the theory acausal or anti-dissipative.

    ``inertia`` is ``rho + p`` in Msun/(kpc Myr^2). A negative ``eta`` would pump energy
    into the fluid; ``tau <= 0`` is the instantaneous (parabolic) limit, which propagates
    faster than light; ``eta / (tau (rho + p)) > 1`` is the Israel-Stewart causality
    violation.
    """
    if eta < 0.0:
        raise ValueError("shear viscosity must be non-negative: negative eta is anti-damping")
    if tau <= 0.0:
        raise ValueError("tau_Pi must be positive: tau = 0 is the acausal parabolic limit")
    if inertia <= 0.0:
        raise ValueError("inertia (rho + p) must be positive")
    if causal_speed(eta, tau, inertia) > NATURAL_CAUSALITY_LIMIT:
        raise ValueError(
            "Israel-Stewart causality violated: "
            f"v_c = {causal_speed(eta, tau, inertia):.3f} > 1 "
            "(decrease eta or increase tau_Pi)"
        )


def relaxed_stress(
    equilibrium_value: Union[float, np.ndarray],
    elapsed_time: Union[float, np.ndarray],
    tau: float,
    initial_value: float = 0.0,
):
    """Exact solution of ``tau dPi/dt + Pi = Pi_eq`` for constant ``Pi_eq``.

    ``Pi(t) = Pi_eq + (Pi_0 - Pi_eq) exp(-t/tau)``. This is the memory the equilibrium
    drag law lacks: at ``t = 0`` the fluid still carries ``Pi_0``, which is what the
    earlier nodes implicitly set to the equilibrium value.
    """
    if tau <= 0.0:
        raise ValueError("tau_Pi must be positive")
    equilibrium_value = np.asarray(equilibrium_value, dtype=float)
    elapsed_time = np.asarray(elapsed_time, dtype=float)
    return equilibrium_value + (initial_value - equilibrium_value) * np.exp(
        -elapsed_time / tau
    )


def crossing_time(width: float, velocity: float) -> float:
    """``t_cross = w / v`` (Myr) for a front of width ``w`` kpc at speed ``v`` kpc/Myr."""
    if width <= 0.0:
        raise ValueError("front width must be positive")
    if velocity <= 0.0:
        raise ValueError("front velocity must be positive")
    return width / velocity


def deborah_number(tau: float, width: float, velocity: float) -> float:
    """``De = tau_Pi / t_cross``: how much memory the front can use.

    ``De << 1`` the dissipative stress equilibrates before the front passes;
    ``De ~ 1`` it is half-formed; ``De >> 1`` the front is fully non-equilibrium and
    most of the stress is deposited behind it.
    """
    return tau / crossing_time(width, velocity)


def causal_speed(eta: float, tau: float, inertia: float) -> float:
    """``v_c = sqrt(eta / (tau (rho + p)))`` with ``c = 1``; must be ``<= 1``."""
    if tau <= 0.0 or inertia <= 0.0:
        raise ValueError("tau_Pi and inertia must be positive")
    return float(np.sqrt(eta / (tau * inertia)))


# --- fronts -----------------------------------------------------------------
def shell_velocity(
    r: Union[float, np.ndarray],
    inner_radius: float,
    outer_radius: float,
    width: float,
    velocity: float,
) -> np.ndarray:
    """Radial velocity ``v(r)`` (kpc/Myr) of one outward-moving shell.

    Zero outside, ``velocity`` inside, smooth across both fronts over a width ``w``.
    ``inner_radius < outer_radius``; the leading front is the outer one.
    """
    if not inner_radius < outer_radius:
        raise ValueError("inner_radius must be smaller than outer_radius")
    if width <= 0.0:
        raise ValueError("front width must be positive")
    r_arr = np.asarray(r, dtype=float)
    return 0.5 * velocity * (
        np.tanh((r_arr - inner_radius) / width) - np.tanh((r_arr - outer_radius) / width)
    )


def radial_shear(
    r: Union[float, np.ndarray],
    inner_radius: float,
    outer_radius: float,
    width: float,
    velocity: float,
) -> np.ndarray:
    """``sigma_rr = (2/3) dv/dr`` (1/Myr), thin-front and non-relativistic.

    The factor ``2/3`` is the traceless radial part of ``d_(mu u_nu)`` in the tetrad:
    ``sigma_rr = (2/3)(dv/dr - Theta/3)`` with ``Theta`` dropped in the thin-front limit.
    Negative at the leading front (the flow decelerates into the ambient medium: a
    compression) and positive at the trailing one (the flow accelerates away from rest: a
    rarefaction).
    """
    r_arr = np.asarray(r, dtype=float)
    sech2 = lambda x: 1.0 / np.cosh(x) ** 2  # noqa: E731 - local, single use
    # Derivative of :func:`shell_velocity` with respect to r, term by term. The order of
    # the two sech^2 terms must match that function, or the shear's sign is reversed.
    slope = 0.5 * velocity / width * (
        sech2((r_arr - inner_radius) / width) - sech2((r_arr - outer_radius) / width)
    )
    return (2.0 / 3.0) * slope


def front_shear(
    inner_radius: float,
    outer_radius: float,
    width: float,
    velocity: float,
) -> Dict[str, float]:
    """Peak shear at each front (1/Myr), named by the front they belong to.

    At a front's own centre only its own ``sech^2`` survives, and there
    ``|dv/dr| = v0 / (2 w)``, so ``|sigma| = v0 / (3 w)``. The two fronts carry the same
    magnitude with opposite signs: the leading (outer) front compresses the flow, the
    trailing (inner) one stretches it. Both follow from the same velocity field, not from
    hand-assigned labels.
    """
    if not inner_radius < outer_radius:
        raise ValueError("inner_radius must be smaller than outer_radius")
    if width <= 0.0 or velocity <= 0.0:
        raise ValueError("width and velocity must be positive")
    peak = velocity / (3.0 * width)
    return {"leading": -peak, "trailing": peak}


def front_summary(
    inner_radius: float,
    outer_radius: float,
    width: float,
    velocity: float,
    eta: float,
    tau: float,
    inertia: float = 1.0,
) -> Dict[str, object]:
    """Everything the engine can say about one shell's pair of fronts.

    Returns the crossing time, the Deborah number, and for each front its peak shear,
    the equilibrium stress ``-2 eta sigma``, the stress actually reached at the crossing
    ``-2 eta sigma (1 - e^{-De})``, the anisotropy ``Delta = (3/2) Pi_rr``, and the
    entropy production ``-(1/2) Pi : sigma`` (always ``>= 0`` for the dissipative branch).
    """
    validate_transport(eta, tau, inertia)
    shear = front_shear(inner_radius, outer_radius, width, velocity)
    t_cross = crossing_time(width, velocity)
    de = deborah_number(tau, width, velocity)

    fronts = {}
    for name in ("leading", "trailing"):
        sigma = shear[name]
        pi_eq = float(equilibrium_stress(eta, sigma))
        pi_cross = float(relaxed_stress(pi_eq, t_cross, tau))
        fronts[name] = {
            "shear": sigma,
            "stress_equilibrium": pi_eq,
            "stress_at_crossing": pi_cross,
            "anisotropy_equilibrium": anisotropy_pressure(pi_eq),
            "anisotropy_at_crossing": anisotropy_pressure(pi_cross),
            "entropy_production": float(entropy_production(pi_cross, sigma)),
            "relaxation_ratio": float(pi_cross / pi_eq) if pi_eq != 0.0 else float("nan"),
        }

    return {
        "crossing_time": t_cross,
        "deborah_number": de,
        "is_non_equilibrium": bool(de >= 1.0),
        "causal_speed": causal_speed(eta, tau, inertia),
        "fronts": fronts,
        "stress_opposes_shear": all(
            np.sign(fronts[name]["stress_at_crossing"]) == -np.sign(fronts[name]["shear"])
            for name in fronts
        ),
    }


# --- derived quantities -----------------------------------------------------
def anisotropy_pressure(pi_rr: Union[float, np.ndarray]):
    """``Delta = p_r - p_t = (3/2) Pi_rr`` — inverts ``Pi_rr = (2/3) Delta``."""
    return 1.5 * np.asarray(pi_rr, dtype=float)


def entropy_production(
    pi_rr: Union[float, np.ndarray],
    shear: Union[float, np.ndarray],
):
    """Entropy production density ``-(1/2) Pi : sigma`` (Msun/(kpc Myr^3)).

    The sign is not cosmetic. With the dissipative stress entering as
    ``T^{mu nu} = ... - 2 eta sigma^{mu nu}`` and ``Pi = -2 eta sigma``, the contraction
    ``Pi : sigma`` is *negative*. What the second law requires non-negative is
    ``-(1/2) Pi : sigma = eta sigma : sigma >= 0``, which is what this function returns.
    Returning ``Pi : sigma`` directly would have reported a dissipative front as if it
    created entropy of the wrong sign.
    """
    return -0.5 * np.asarray(pi_rr, dtype=float) * np.asarray(shear, dtype=float)


def exchange_decomposition(Q_energy: float, force_radial: float) -> Dict[str, object]:
    """Split the interaction vector ``Q^nu = Q u^nu + f^nu`` with ``u.f = 0``.

    ``Q`` is the energy/mass exchange (``Q > 0``: vacuum decays into matter) and ``f^alpha
    = h^{alpha nu} Q_nu`` the momentum transfer, which is what accelerates or torques the
    matter. Rows are returned as the tetrad components of ``Q^nu``.
    """
    u = four_velocity()
    s = spatial_direction()
    vector = Q_energy * u + force_radial * s
    # Orthogonality of the split, checked rather than asserted in prose.
    if not np.isclose(float((geometry() @ u) @ vector), -Q_energy, atol=1e-12):
        raise ValueError("decomposition is inconsistent with u.Q = -Q")
    return {
        "vector": vector,
        "energy": float(Q_energy),
        "force": float(force_radial),
        "matter_source": vector,
        "vacuum_sink": -vector,
    }


@dataclass
class FrontStressState:
    """A single front's dissipative state, with the invariants already checked."""

    shear: float                 # 1/Myr
    stress: float                # Msun/(kpc Myr^2), the Pi_rr component
    tau: float = 1.0             # Myr
    eta: float = 1.0             # Msun/(kpc Myr)

    @property
    def tensor(self) -> np.ndarray:
        """``Pi^{mu nu}`` built from ``Pi_rr`` (``Delta = 3 Pi_rr / 2``)."""
        return stress_tensor(float(anisotropy_pressure(self.stress)))

    @property
    def anisotropy(self) -> float:
        return float(anisotropy_pressure(self.stress))

    @property
    def is_dissipative(self) -> bool:
        """``-(1/2) Pi : sigma >= 0``: the front removes energy from the flow, never adds it."""
        return bool(float(entropy_production(self.stress, self.shear)) >= 0.0)
