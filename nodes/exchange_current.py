"""
nodes/exchange_current.py
The interaction four-vector ``J^mu`` and the total conservation law that fixes it.

Why this module exists
----------------------
Node 05 writes the exchange as ``J^mu = Gamma (u_v^mu - u_m^mu)`` and then keeps two
numbers from it: a drag force ``-Gamma du`` and an energy rate ``-Gamma du^2``. In that
form the object is *radial and scalar*: four components in principle, one in practice, and
the azimuthal one — the one that carries torque — is dropped before it is ever computed.

This module keeps the four. For an arbitrary frame it builds

    J^mu = Q u^mu + f^mu,          u_mu f^mu = 0

with ``Q = -u_mu J^mu`` the energy exchange and ``f^alpha = h^alpha{}_nu J^nu`` the momentum
transfer. The split is not a convention applied afterwards: :class:`ExchangeCurrent`
refuses a force with a component along ``u``, and :meth:`ExchangeCurrent.from_vector`
recovers the split of *any* four-vector exactly, so nothing is lost on the way in.

Conservation
------------
``nabla_mu T_tot^{mu nu} = 0`` is a Bianchi identity, not an assumption, so the matter
source is fixed by the vacuum sector rather than chosen:

    nabla_mu T_m^{mu nu} = -nabla_mu T_vac^{mu nu} = J^nu.

For the anisotropic vacuum stress

    T_vac^{mu nu} = (rho_v + p_t) u^mu u^nu + p_t g^{mu nu} + Delta s^mu s^nu + Pi^{mu nu},

contracting the divergence with ``u`` is an identity, not a model statement:

    u_nu nabla_mu T_vac^{mu nu} = -[ rho_dot_v + (rho_v + p_t) Theta
                                     + Delta (s^mu s^nu nabla_mu u_nu) + Pi : sigma ]

so

    Q = rho_dot_v + (rho_v + p_t) Theta + Delta (s:grad u:s) + Pi:sigma.  (*)

``Theta = nabla.u`` is the expansion, ``Pi:sigma`` the dissipative contraction. Three
things in (*) are worth naming. The pressure that appears is ``p_t``, not the isotropic
mean: an anisotropic vacuum converts differently along the two axes. The ``Delta`` term is
the *stretching of the anisotropy direction along itself*, which a purely radial model sets
to zero by construction. And ``Pi:sigma`` remains even where ``Delta = 0``, because the
dissipative stress is its own source.

:func:`vacuum_energy_exchange` returns ``Q`` from (*), and the suite checks it against a
finite-difference divergence of the same stress tensor. That is a real test: drop the
``Pi`` term, or replace ``p_t`` with the mean pressure, and the two sides separate at
second order in the grid spacing.

Flat-space grid helpers
-----------------------
:func:`vacuum_stress`, :func:`flat_divergence` and :func:`conservation_ledger` operate on
a ``(t, r)`` grid in an orthonormal frame with signature ``(-,+,+,+)`` and
``x^mu = (t, r, theta, phi)``. Fluxes in ``theta`` and ``phi`` are zero by axisymmetry;
the azimuthal *velocity* is not, and it is what makes ``f^phi`` non-zero. No metric factors
enter because the basis is orthonormal, which is also why no Christoffel symbols appear.

Units
-----
length kpc, time Myr, velocity kpc/Myr, stress and pressure Msun/(kpc Myr^2), four-current
``J^mu`` Msun/(kpc^2 Myr^2) — a stress per length, i.e. a force density — torque density
``r f^phi`` Msun/(kpc Myr^2). ``c = 1`` throughout.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Union

import numpy as np

from .non_equilibrium_stress import four_velocity, geometry, shell_velocity

__all__ = [
    "validate_frame",
    "four_velocity",
    "boosted_velocity",
    "radial_direction",
    "spatial_projector",
    "flow_projector",
    "project_onto_flow_orthogonal",
    "ExchangeCurrent",
    "exchange_from_relative_velocity",
    "vacuum_stress",
    "vacuum_stress_energy_in_frame",
    "flat_divergence",
    "vacuum_force",
    "euler_acceleration",
    "vacuum_energy_exchange",
    "torque_density",
    "conservation_ledger",
    "rotating_front_current",
]

#: Absolute tolerance for the frame conditions (unit vectors, orthogonality).
FRAME_ATOL = 1e-10


# --- inner products ---------------------------------------------------------
def _dot(a, b):
    """``a.b`` with signature ``(-,+,+,+)``, broadcasting over leading axes."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return -a[..., 0] * b[..., 0] + np.sum(a[..., 1:] * b[..., 1:], axis=-1)


def _outer(a, b):
    """``a^mu b^nu`` with broadcasting over leading axes."""
    return np.einsum("...i,...j->...ij", np.asarray(a, dtype=float), np.asarray(b, dtype=float))


# --- frame ------------------------------------------------------------------
def validate_frame(u=None, s=None) -> tuple:
    """Enforce ``u.u = -1``, ``s.s = +1``, ``u.s = 0``; return the pair.

    Accepts arrays of shape ``(..., 4)`` so a whole grid is checked at once. A vector that
    is merely *normalized* is not enough: ``s`` must also be orthogonal to ``u``, otherwise
    ``s^mu`` is not a rest-frame axis and ``T^{mu nu}`` carries a rest-frame energy flux.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    if np.any(~np.isfinite(u)):
        raise ValueError("four-velocity must be finite")
    norm_u = _dot(u, u)
    if np.any(np.abs(norm_u + 1.0) > FRAME_ATOL):
        raise ValueError(
            f"u.u = {np.max(np.abs(norm_u)):.6g} off, expected -1 (timelike, (-,+,+,+))"
        )
    if s is None:
        return u, None
    s = np.asarray(s, dtype=float)
    if np.any(~np.isfinite(s)):
        raise ValueError("spatial direction must be finite")
    norm_s = _dot(s, s)
    if np.any(np.abs(norm_s - 1.0) > FRAME_ATOL):
        raise ValueError(
            f"s.s = {np.min(np.abs(norm_s)):.6g} off, expected +1 (spacelike unit vector)"
        )
    overlap = _dot(u, s)
    if np.any(np.abs(overlap) > FRAME_ATOL):
        raise ValueError(
            f"u.s = {np.max(np.abs(overlap)):.6g}: the radial direction must be "
            "orthogonal to the flow"
        )
    return u, s


def boosted_velocity(beta_r: float = 0.0, beta_theta: float = 0.0, beta_phi: float = 0.0):
    """``u^mu = gamma (1, beta_r, beta_theta, beta_phi)`` in the orthonormal tetrad.

    ``|beta| < 1`` is enforced: a null or superluminal flow has no rest frame, so ``u`` is
    not a four-velocity and every contraction below would be meaningless.
    """
    beta = np.array([beta_r, beta_theta, beta_phi], dtype=float)
    if np.any(~np.isfinite(beta)):
        raise ValueError("velocity components must be finite")
    beta_sq = float(beta @ beta)
    if beta_sq >= 1.0:
        raise ValueError(f"|beta|^2 = {beta_sq:.6g} >= 1: subluminal flow required")
    gamma = 1.0 / np.sqrt(1.0 - beta_sq)
    return np.concatenate(([gamma], gamma * beta))


def radial_direction(u=None):
    """Unit radial axis of the *rest frame* of ``u``: ``s = (e_r + (u.e_r) u)/sqrt(1+u_r^2)``.

    For a flow at rest this is ``e_r``. For a radially moving flow it is unchanged, because
    ``e_r`` was already orthogonal to ``u``. For a flow with an azimuthal component it picks
    up a time component, and that component is not cosmetic: it is what keeps
    ``T^{0i} = 0`` in the rest frame. Handing ``e_r`` to a rotating flow would leave a
    rest-frame energy flux in an equilibrium configuration.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    u, _ = validate_frame(u)
    e_r = np.zeros(np.shape(u))
    e_r[..., 1] = 1.0
    u_r = _dot(u, e_r)
    raw = e_r + np.expand_dims(u_r, -1) * u
    norm = np.sqrt(1.0 + u_r**2)
    return raw / np.expand_dims(norm, -1)


def spatial_projector(u=None):
    """``h^{mu nu} = g^{mu nu} + u^mu u^nu`` for an arbitrary ``u`` (not only ``(1,0,0,0)``).

    This is the *all-up* projector. It acts on a **covector**: ``h^{mu nu} w_nu``. Applying
    it to a contravariant vector with a matrix multiplication is a different object and
    gives a non-zero result for ``w = u``; use :func:`project_onto_flow_orthogonal` for
    that, which is the mixed projector ``h^mu{}_nu = delta^mu_nu + u^mu u_nu``.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    u, _ = validate_frame(u)
    return geometry() + _outer(u, u)


def flow_projector(u=None):
    """Mixed projector ``h^mu{}_nu = delta^mu_nu + u^mu u_nu`` — the one that acts on vectors."""
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    u, _ = validate_frame(u)
    lowered = np.einsum("...j,jk->...k", u, geometry())
    return np.eye(4) + _outer(u, lowered)


def project_onto_flow_orthogonal(vector, u=None):
    """``h^alpha{}_nu v^nu = v^alpha + u^alpha (u.v)``: the part of ``v`` orthogonal to ``u``.

    Uses ``u.v`` (with the metric), not an unweighted dot: writing ``h^{alpha nu} v_nu``
    and then multiplying by a contravariant ``v`` is the classic place where a projection
    that should vanish at ``v = u`` does not.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    u, _ = validate_frame(u)
    vector = np.asarray(vector, dtype=float)
    if vector.shape != u.shape:
        raise ValueError(f"shape mismatch: {vector.shape} vs {u.shape}")
    return vector + np.expand_dims(_dot(u, vector), -1) * u


# --- the current ------------------------------------------------------------
@dataclass
class ExchangeCurrent:
    """``J^mu = Q u^mu + f^mu`` with ``u.f = 0``.

    ``Q`` is the energy exchange: positive means the vacuum sector loses energy to the
    matter sector. ``f`` is the momentum transfer, which is what accelerates *and torques*
    the matter. The force is stored contravariant, so ``f^phi`` is a component in the
    orthonormal frame and the torque density is ``r f^phi``.
    """

    energy: float
    force: np.ndarray
    u: np.ndarray = field(default_factory=four_velocity)

    def __post_init__(self):
        self.u = np.asarray(self.u, dtype=float).reshape(4)
        self.force = np.asarray(self.force, dtype=float).reshape(4)
        self.energy = float(self.energy)
        if not np.isfinite(self.energy) or np.any(~np.isfinite(self.force)):
            raise ValueError("exchange current must be finite")
        validate_frame(self.u)
        overlap = float(_dot(self.force, self.u))
        if abs(overlap) > FRAME_ATOL:
            raise ValueError(
                f"f.u = {overlap:.6g}: the momentum transfer must be orthogonal to the "
                "flow; it is h^{alpha nu} J_nu, not an unconverted four-vector"
            )

    # -- construction --------------------------------------------------------
    @classmethod
    def from_vector(cls, vector, u=None) -> "ExchangeCurrent":
        """Split an arbitrary ``J^nu`` into ``Q = -u.J`` and ``f = h J``; returns the same vector."""
        u = four_velocity() if u is None else np.asarray(u, dtype=float).reshape(4)
        vector = np.asarray(vector, dtype=float).reshape(4)
        validate_frame(u)
        energy = -float(_dot(u, vector))
        return cls(energy=energy, force=project_onto_flow_orthogonal(vector, u), u=u)

    # -- components ----------------------------------------------------------
    @property
    def vector(self) -> np.ndarray:
        """``J^mu`` in the frame of ``u``."""
        return self.energy * self.u + self.force

    @property
    def matter_source(self) -> np.ndarray:
        """What the matter sector gains: ``+J^mu``."""
        return self.vector

    @property
    def vacuum_sink(self) -> np.ndarray:
        """What the vacuum sector loses: ``-J^mu``. The pair sums to zero exactly."""
        return -self.vector

    @property
    def momentum(self) -> np.ndarray:
        """``h^alpha{}_nu J^nu = f^alpha``, recovered from the vector rather than stored."""
        return project_onto_flow_orthogonal(self.vector, self.u)

    @property
    def norm(self) -> float:
        """``|J|^2 = f.f - Q^2`` (with ``c = 1``)."""
        return float(_dot(self.vector, self.vector))

    @property
    def is_pure_energy(self) -> bool:
        return bool(np.all(np.abs(self.force) <= FRAME_ATOL))

    @property
    def is_pure_force(self) -> bool:
        return abs(self.energy) <= FRAME_ATOL

    def components(self) -> Dict[str, float]:
        return {
            "energy": self.energy,
            "f_r": float(self.force[1]),
            "f_theta": float(self.force[2]),
            "f_phi": float(self.force[3]),
        }

    # -- torque --------------------------------------------------------------
    def torque_density(self, radius: float) -> float:
        """``r f^phi``: angular-momentum transfer per volume per time (Msun/(kpc Myr^2)).

        A purely radial force has ``f^phi = 0`` and exerts no torque, which is exactly the
        sector the radial-only implementation could not reach.
        """
        if radius < 0.0:
            raise ValueError("radius must be non-negative")
        return float(radius) * float(self.force[3])

    def angular_momentum_rate(self, radius: float) -> float:
        """Alias of :meth:`torque_density`, named for what it conserves."""
        return self.torque_density(radius)


def exchange_from_relative_velocity(gamma_coefficient, u_matter, u_vacuum) -> ExchangeCurrent:
    """``J^mu = Gamma (u_v^mu - u_m^mu)`` — Node 05's law, with its full components.

    Node 05 keeps only ``F = -Gamma du`` and ``u_m.J = -Gamma (du)^2``, which are the
    ``O(beta)`` and ``O(beta^2)`` pieces of the same current. Here the whole vector is
    available, so an azimuthal drift between the two fluids produces ``f^phi`` and
    therefore a torque, without any extra assumption.
    """
    gamma_coefficient = float(gamma_coefficient)
    if not np.isfinite(gamma_coefficient):
        raise ValueError("exchange coefficient must be finite")
    if gamma_coefficient < 0.0:
        raise ValueError("exchange coefficient must be non-negative: it cannot act as a drive")
    u_matter = np.asarray(u_matter, dtype=float).reshape(4)
    u_vacuum = np.asarray(u_vacuum, dtype=float).reshape(4)
    validate_frame(u_matter)
    validate_frame(u_vacuum)
    return ExchangeCurrent.from_vector(gamma_coefficient * (u_vacuum - u_matter), u=u_matter)


# --- the vacuum stress tensor ----------------------------------------------
def vacuum_stress(
    rho: Union[float, np.ndarray],
    p_t: Union[float, np.ndarray],
    delta: Union[float, np.ndarray] = 0.0,
    u=None,
    s=None,
    pi_amplitude: Union[float, np.ndarray] = 0.0,
):
    """``T_vac^{mu nu}`` for an anisotropic, dissipative vacuum.

    ``T = (rho + p_t) u u + p_t g + Delta s s + Pi`` with ``Pi = pi_amplitude * chi``,
    i.e. ``Pi_rr = (2/3) pi_amplitude``. Accepts scalars or grids; ``u`` and ``s`` may have
    a leading grid axis. ``s`` defaults to the rest-frame radial axis of ``u``.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    if np.ndim(u) == 1:
        u = u.reshape(4)
        if s is None:
            s = radial_direction(u)
        s = np.asarray(s, dtype=float).reshape(4)
        validate_frame(u, s)
        h = geometry() + np.outer(u, u)
        tensor = (
            (float(rho) + float(p_t)) * np.outer(u, u)
            + float(p_t) * geometry()
            + float(np.asarray(delta)) * np.outer(s, s)
            + float(np.asarray(pi_amplitude)) * (np.outer(s, s) - h / 3.0)
        )
        return tensor

    u = u.reshape(u.shape)
    if s is None:
        s = radial_direction(u)
    s = np.asarray(s, dtype=float)
    validate_frame(u, s)
    h = geometry() + _outer(u, u)
    scale = lambda value: np.asarray(value, dtype=float)[..., None, None]  # noqa: E731
    return (
        scale(np.asarray(rho) + np.asarray(p_t)) * _outer(u, u)
        + scale(p_t) * geometry()
        + scale(delta) * _outer(s, s)
        + scale(pi_amplitude) * (_outer(s, s) - h / 3.0)
    )


def vacuum_stress_energy_in_frame(rho, p_t, delta=0.0, u=None, s=None, pi_amplitude=0.0):
    """``u_nu T^{mu nu} = -rho u^mu``: the anisotropic stress carries no rest-frame energy.

    This is exact for the tensor above, and it is the sharpest reason the anisotropy
    direction must be orthogonal to ``u``. If ``s`` is merely normalized, ``u.s != 0`` and
    the right-hand side stops being ``-rho u^mu``.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float).reshape(4)
    tensor = vacuum_stress(rho, p_t, delta=delta, u=u, s=s, pi_amplitude=pi_amplitude)
    return (geometry() @ u) @ tensor


# --- derivatives on a flat grid --------------------------------------------
def flat_divergence(field, dt: float, dr: float):
    """``d_mu T^{mu nu}`` on a ``(t, r)`` grid in the orthonormal frame.

    Only ``mu = t`` and ``mu = r`` contribute: the configuration is axisymmetric, so
    ``d_theta`` and ``d_phi`` vanish. ``field`` has shape ``(nt, nr, 4, 4)``; the result has
    shape ``(nt, nr, 4)``. Second-order central differences interior, one-sided at the
    boundary, so interior points are the accurate ones.
    """
    field = np.asarray(field, dtype=float)
    if field.ndim != 4 or field.shape[2:] != (4, 4):
        raise ValueError("field must have shape (nt, nr, 4, 4)")
    if dt <= 0.0 or dr <= 0.0:
        raise ValueError("grid spacing must be positive")
    d_time = np.gradient(field, dt, axis=0)
    d_radius = np.gradient(field, dr, axis=1)
    return d_time[:, :, 0, :] + d_radius[:, :, 1, :]


def vacuum_force(divergence, u=None):
    """``F_vac^alpha = -h^alpha{}_nu d_mu T_vac^{mu nu}``: the force on the matter.

    Isotropic and uniform vacuum gives ``d T = 0`` and no local force, whatever the energy
    density. With ``Pi != 0`` the bracket survives even where ``Delta = 0``, because the
    dissipative stress is itself a gradient source.
    """
    u = four_velocity() if u is None else np.asarray(u, dtype=float)
    divergence = np.asarray(divergence, dtype=float)
    validate_frame(u)
    return -project_onto_flow_orthogonal(divergence, u)


def euler_acceleration(inertia, pressure_gradient, force, u=None):
    """``a^alpha = (f^alpha - h^{alpha mu} d_mu p_m) / (rho_m + p_m)``.

    The generalized Euler equation with the vacuum as the only source. ``inertia`` is
    ``rho_m + p_m`` (with ``c = 1``), ``pressure_gradient`` is ``d_mu p_m`` in the same
    orthonormal frame, ``force`` is ``f^alpha`` from :class:`ExchangeCurrent`. The result is
    orthogonal to ``u`` by construction, because both terms are projected.
    """
    inertia = float(inertia)
    if not np.isfinite(inertia) or inertia <= 0.0:
        raise ValueError("inertia (rho_m + p_m) must be positive and finite")
    u = four_velocity() if u is None else np.asarray(u, dtype=float).reshape(4)
    validate_frame(u)
    gradient = np.asarray(pressure_gradient, dtype=float).reshape(4)
    force = np.asarray(force, dtype=float).reshape(4)
    if abs(float(_dot(force, u))) > FRAME_ATOL:
        raise ValueError("the vacuum force must be orthogonal to the flow")
    return (force - project_onto_flow_orthogonal(gradient, u)) / inertia


# --- energy balance ---------------------------------------------------------
def vacuum_energy_exchange(
    rho_dot: float,
    rho: float,
    p_t: float,
    theta: float,
    delta: float = 0.0,
    direction_stretch: float = 0.0,
    pi_sigma: float = 0.0,
) -> float:
    """``Q = rho_dot + (rho + p_t) Theta + Delta (s:grad u:s) + Pi:sigma``.

    The identity obtained by contracting ``nabla_mu T_vac^{mu nu} = -J^nu`` with ``u_nu``.
    ``Q > 0``: the vacuum sector loses energy to the matter sector.

    Each term is separately meaningful. ``rho_dot`` is the rest-frame density change;
    ``(rho + p_t) Theta`` is the ``p dV`` work, with the *tangential* pressure because that
    is the one the isotropic part of the anisotropy leaves behind; ``direction_stretch``
    is how much the anisotropy axis is stretched along itself; and ``pi_sigma`` is the
    dissipative contraction, which survives even where ``Delta = 0``.
    """
    for name, value in (
        ("rho_dot", rho_dot),
        ("rho", rho),
        ("p_t", p_t),
        ("theta", theta),
        ("delta", delta),
        ("direction_stretch", direction_stretch),
        ("pi_sigma", pi_sigma),
    ):
        if not np.isfinite(float(value)):
            raise ValueError(f"{name} must be finite")
    return float(
        rho_dot
        + (rho + p_t) * theta
        + delta * direction_stretch
        + pi_sigma
    )


def torque_density(radius: float, f_phi: float) -> float:
    """``r f^phi`` in Msun/(kpc Myr^2), the same quantity as :meth:`ExchangeCurrent.torque_density`."""
    radius = float(radius)
    if radius < 0.0:
        raise ValueError("radius must be non-negative")
    if not np.isfinite(radius) or not np.isfinite(float(f_phi)):
        raise ValueError("inputs must be finite")
    return radius * float(f_phi)


# --- conservation ledger ----------------------------------------------------
def conservation_ledger(matter_divergence, vacuum_divergence) -> Dict[str, object]:
    """Check ``nabla.T_m + nabla.T_vac = 0`` and report the residual.

    The two sectors are supplied independently — normally one is measured from the matter
    and the other from the vacuum stress — so a sign error in the exchange convention shows
    up as a non-zero residual instead of cancelling silently. ``residual`` is the norm of
    the total divergence; ``balanced`` is ``residual <= atol``.
    """
    matter = np.asarray(matter_divergence, dtype=float)
    vacuum = np.asarray(vacuum_divergence, dtype=float)
    if matter.shape != vacuum.shape:
        raise ValueError(f"shape mismatch: {matter.shape} vs {vacuum.shape}")
    if matter.shape[-1] != 4:
        raise ValueError("a divergence of T^{mu nu} has four components")
    total = matter + vacuum
    residual = float(np.max(np.abs(total))) if total.size else 0.0
    atol = 1e-12 * max(1.0, float(np.max(np.abs(matter))) + float(np.max(np.abs(vacuum))))
    return {
        "total": total,
        "residual": residual,
        "balanced": bool(residual <= atol),
        "matter_magnitude": float(np.max(np.abs(matter))),
        "vacuum_magnitude": float(np.max(np.abs(vacuum))),
    }


# --- a worked front, for the documentation table ----------------------------
def rotating_front_current(
    r: float,
    inner_radius: float,
    outer_radius: float,
    width: float,
    radial_velocity: float,
    azimuthal_velocity: float,
    gamma_coefficient: float,
):
    """``J^mu`` across a front where the vacuum lags the matter azimuthally.

    The vacuum is at rest, the matter carries ``v_r(r)`` from :func:`shell_velocity` and a
    constant azimuthal ``v_phi``. Returns the current and the torque density at ``r``.
    Used only to show that the azimuthal sector exists and has the right sign; no parameter
    is fitted.
    """
    v_r = float(shell_velocity(r, inner_radius, outer_radius, width, radial_velocity))
    beta = np.array([v_r, 0.0, azimuthal_velocity], dtype=float)
    if float(beta @ beta) >= 1.0:
        raise ValueError("|beta| >= 1: subluminal flow required")
    u_matter = boosted_velocity(*beta)
    current = exchange_from_relative_velocity(gamma_coefficient, u_matter, four_velocity())
    return current, current.torque_density(r)

