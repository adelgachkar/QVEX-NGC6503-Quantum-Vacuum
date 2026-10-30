"""Invariants of the interaction four-vector ``J^mu`` and of total conservation.

The subject of these tests is the difference between a *radial and scalar* exchange law —
``F = -Gamma du`` with one number at each radius, as Node 05 keeps it — and the full
four-vector ``J^mu = Q u^mu + f^mu``. Two tests carry that argument:

* ``test_decomposition_keeps_the_azimuthal_component`` — a current with a ``phi`` component
  survives the split, so the torque sector is reachable at all;
* ``test_energy_identity_matches_the_numerical_divergence`` — the identity
  ``Q = rho_dot + (rho + p_t) Theta + Delta (s:grad u:s) + Pi:sigma`` is checked against a
  finite-difference divergence of the same stress tensor. Drop the ``Pi`` term, or use the
  mean pressure instead of ``p_t``, and the two sides separate at second order in ``h``.

No golden numbers: every assertion is a property that must hold for any parameters.
"""

import numpy as np
import pytest

from nodes.exchange_current import (
    FRAME_ATOL,
    ExchangeCurrent,
    _dot,
    _outer,
    boosted_velocity,
    flow_projector,
    project_onto_flow_orthogonal,
    conservation_ledger,
    euler_acceleration,
    exchange_from_relative_velocity,
    flat_divergence,
    radial_direction,
    spatial_projector,
    torque_density,
    vacuum_energy_exchange,
    vacuum_force,
    vacuum_stress,
    vacuum_stress_energy_in_frame,
    validate_frame,
)
from nodes.non_equilibrium_stress import four_velocity, geometry

METRIC = geometry()


# --- frame conditions -------------------------------------------------------
@pytest.mark.parametrize("beta_r", [0.0, 0.37, -0.6])
@pytest.mark.parametrize("beta_phi", [0.0, 0.42])
def test_rest_frame_axes_are_orthonormal_for_a_boosted_and_rotating_flow(beta_r, beta_phi):
    """u.u = -1, s.s = +1, u.s = 0 — for an oblique flow, not only for (1,0,0,0)."""
    u = boosted_velocity(beta_r, 0.0, beta_phi)
    s = radial_direction(u)

    assert np.isclose(float(_dot(u, u)), -1.0)
    assert np.isclose(float(_dot(s, s)), 1.0)
    assert np.isclose(float(_dot(u, s)), 0.0, atol=1e-13)
    # A radially moving flow tilts the radial axis into the time direction; that tilt is what
    # keeps the rest-frame energy flux zero, so it must not be silently dropped.
    if abs(beta_r) > 0.0:
        assert abs(s[0]) > 0.0
    else:
        assert s[0] == 0.0      # the flow at rest, or purely rotating, is orthogonal to e_r


def test_superluminal_and_unnormalized_frames_are_rejected():
    with pytest.raises(ValueError, match="subluminal"):
        boosted_velocity(1.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="subluminal"):
        boosted_velocity(0.8, 0.0, 0.8)

    with pytest.raises(ValueError, match="timelike"):
        validate_frame(np.array([1.0, 0.1, 0.0, 0.0]))

    # Normalized but not orthogonal: this is the failure a check on |s| alone would miss.
    tilted = radial_direction(boosted_velocity(0.3, 0.0, 0.4)) + 0.5 * four_velocity()
    tilted = tilted / np.sqrt(float(_dot(tilted, tilted)))
    assert np.isclose(float(_dot(tilted, tilted)), 1.0)
    with pytest.raises(ValueError, match="orthogonal to the flow"):
        validate_frame(four_velocity(), tilted)


# --- the current ------------------------------------------------------------
def test_decomposition_keeps_the_azimuthal_component():
    """A radial-only implementation loses ``f^phi``; the split must not."""
    u = boosted_velocity(0.2, 0.0, 0.3)
    h = spatial_projector(u)
    tangential = np.array([0.0, 0.0, 0.0, 1.0])
    vector = 0.5 * u + tangential          # arbitrary, but with a phi component
    current = ExchangeCurrent.from_vector(vector, u)

    assert np.allclose(current.vector, vector, atol=1e-12)   # nothing lost on the way in
    assert np.isclose(float(_dot(current.force, u)), 0.0, atol=1e-13)
    assert np.allclose(current.momentum, project_onto_flow_orthogonal(vector, u), atol=1e-12)
    # e_phi is not itself orthogonal to a rotating flow, but the projected force keeps a
    # large phi part — that is the component a radial model would report as zero.
    assert abs(current.force[3]) > 1.0
    assert abs(current.force[1]) < 0.1

    # And it is not a vector along the flow: a radial model cannot produce this.
    assert not np.isclose(current.torque_density(10.0), 0.0)


def test_a_force_along_the_flow_is_refused():
    """``u.f = 0`` is a definition of ``f``, not a numerical coincidence."""
    u = four_velocity()
    with pytest.raises(ValueError, match="orthogonal to the flow"):
        ExchangeCurrent(energy=1.0, force=np.array([0.3, 1.0, 0.0, 0.0]), u=u)
    with pytest.raises(ValueError, match="finite"):
        ExchangeCurrent(energy=np.inf, force=np.zeros(4))


def test_energy_and_force_sectors_are_distinguishable():
    u = four_velocity()
    pure_energy = ExchangeCurrent(energy=2.0, force=np.zeros(4), u=u)
    pure_force = ExchangeCurrent(energy=0.0, force=np.array([0.0, 1.0, 0.0, 0.0]), u=u)

    assert pure_energy.is_pure_energy and not pure_energy.is_pure_force
    assert pure_force.is_pure_force and not pure_energy.is_pure_force
    assert np.isclose(float(_dot(pure_energy.vector, u)), -2.0)
    assert pure_energy.components() == {"energy": 2.0, "f_r": 0.0, "f_theta": 0.0, "f_phi": 0.0}
    assert pure_force.components()["f_r"] == 1.0
    assert np.isclose(pure_force.norm, 1.0)     # f.f - Q^2 with c = 1


def test_matter_and_vacuum_legs_sum_to_zero_exactly():
    u = boosted_velocity(0.1, 0.0, 0.2)
    force = project_onto_flow_orthogonal(np.array([0.0, 0.4, 0.0, 0.2]), u)
    current = ExchangeCurrent(energy=0.7, force=force, u=u)

    assert np.allclose(current.matter_source + current.vacuum_sink, 0.0)
    ledger = conservation_ledger(current.matter_source, current.vacuum_sink)
    assert ledger["balanced"] is True and ledger["residual"] == 0.0

    # The check is real: the same current on both legs does not balance.
    assert conservation_ledger(current.vector, current.vector)["balanced"] is False


# --- the reduction of Node 05 ----------------------------------------------
def test_relative_velocity_current_reproduces_the_scalar_law_at_low_speed():
    """Node 05's ``F = -Gamma du`` and ``u.J = -Gamma du^2`` are the O(beta) and O(beta^2) parts."""
    gamma_coefficient = 1.7
    dv = 1e-3
    u_matter = boosted_velocity(dv, 0.0, 0.0)
    current = exchange_from_relative_velocity(gamma_coefficient, u_matter, four_velocity())

    # Energy: u_m.J = Gamma (1 - gamma) < 0, so Q = -u_m.J = Gamma dv^2 / 2 + O(dv^4).
    energy_rate = float(_dot(u_matter, gamma_coefficient * (four_velocity() - u_matter)))
    assert energy_rate < 0.0                        # matter gives up energy to the vacuum...
    assert current.energy > 0.0                     # ...and Q = -u.J records it as a gain
    assert np.isclose(current.energy, gamma_coefficient * (1.0 / np.sqrt(1 - dv**2) - 1.0),
                      rtol=1e-12)
    assert np.isclose(current.energy, gamma_coefficient * dv**2 / 2.0, rtol=1e-6)
    assert np.isclose(current.energy, -energy_rate, rtol=1e-12)

    # Force: f^r = -Gamma gamma dv, i.e. -Gamma du in the non-relativistic limit: a drag.
    assert np.isclose(current.force[1], -gamma_coefficient * dv / np.sqrt(1 - dv**2))

    with pytest.raises(ValueError, match="non-negative"):
        exchange_from_relative_velocity(-1.0, u_matter, four_velocity())


def test_azimuthal_drift_becomes_a_torque_with_the_right_sign():
    """A vacuum that lags a rotating flow brakes it; this is the sector Node 05 cannot express."""
    u_matter = boosted_velocity(0.0, 0.0, 0.4)
    lagging = exchange_from_relative_velocity(1.0, u_matter, four_velocity())
    leading = exchange_from_relative_velocity(1.0, four_velocity(), u_matter)

    assert lagging.torque_density(8.0) < 0.0     # f^phi < 0: the matter is slowed
    assert leading.torque_density(8.0) > 0.0     # vacuum ahead: the matter is spun up
    # The two cases are the same physical current seen from two rest frames, so the torques
    # differ by the boost factor gamma of the rotating flow -- not by an arbitrary constant.
    gamma = 1.0 / np.sqrt(1.0 - 0.4**2)
    assert np.isclose(lagging.torque_density(8.0), -gamma * leading.torque_density(8.0))

    # Both currents are orthogonal to their own flow: the torque is a projection, not a leak.
    for matter, current in ((u_matter, lagging), (four_velocity(), leading)):
        assert np.isclose(float(_dot(current.force, matter)), 0.0, atol=1e-13)
    # A purely radial force exerts no torque at all.
    radial = exchange_from_relative_velocity(1.0, boosted_velocity(0.3), four_velocity())
    assert radial.torque_density(12.0) == 0.0
    assert np.isclose(torque_density(12.0, 0.25), 3.0)


# --- the stress tensor -----------------------------------------------------
@pytest.mark.parametrize("delta", [-1.0, 0.0, 2.5])
@pytest.mark.parametrize("pi_amplitude", [0.0, 0.7])
def test_anisotropy_carries_no_energy_in_the_rest_frame(delta, pi_amplitude):
    """``u.T = -rho u`` exactly: the anisotropy moves momentum, never rest-frame energy."""
    u = boosted_velocity(0.25, 0.0, 0.35)
    s = radial_direction(u)
    rho, p_t = 1.3, 0.4

    contracted = vacuum_stress_energy_in_frame(rho, p_t, delta, u=u, s=s,
                                               pi_amplitude=pi_amplitude)
    assert np.allclose(contracted, -rho * u, atol=1e-12)


def test_the_projector_and_the_force_are_orthogonal_to_the_flow():
    u = boosted_velocity(0.3, 0.0, 0.2)
    h_up = spatial_projector(u)
    h_mixed = flow_projector(u)

    # h^mu{}_nu u^nu = 0 is the statement that matters; h^{mu nu} u_nu = 0 is its covector twin.
    assert np.allclose(h_mixed @ u, 0.0, atol=1e-13)
    assert np.allclose(h_up @ (METRIC @ u), 0.0, atol=1e-13)
    # h^{mu nu} is traced with the metric (it is all-up), h^mu{}_nu without it.
    assert np.isclose(float(np.trace(METRIC @ h_up)), 3.0)
    assert np.isclose(float(np.trace(h_mixed)), 3.0)

    arbitrary = np.array([0.4, -0.7, 0.2, 1.1])
    force = vacuum_force(arbitrary, u)
    assert np.isclose(float(_dot(force, u)), 0.0, atol=1e-13)
    # Uniform isotropic vacuum: the divergence vanishes, so no local force.
    assert np.allclose(vacuum_force(np.zeros(4), u), 0.0)


def test_euler_acceleration_reduces_to_the_pressure_gradient_and_to_the_force():
    u = boosted_velocity(0.2, 0.0, 0.1)
    gradient = np.array([0.0, 0.3, 0.0, 0.0])
    inertia = 2.0

    pressure_only = euler_acceleration(inertia, gradient, np.zeros(4), u)
    assert np.isclose(float(_dot(pressure_only, u)), 0.0, atol=1e-13)
    assert pressure_only[1] < 0.0            # pressure pushes outward, the gradient pulls in
    assert np.allclose(pressure_only, -project_onto_flow_orthogonal(gradient, u) / inertia)

    force = project_onto_flow_orthogonal(np.array([0.0, 0.0, 0.0, 0.5]), u)
    pushed = euler_acceleration(inertia, np.zeros(4), force, u)
    assert pushed[3] > 0.0                   # tangential force -> tangential acceleration
    assert np.allclose(pushed, force / inertia)
    assert np.isclose(float(_dot(pushed, u)), 0.0, atol=1e-13)

    with pytest.raises(ValueError, match="positive"):
        euler_acceleration(0.0, gradient, np.zeros(4), u)
    with pytest.raises(ValueError, match="orthogonal"):
        euler_acceleration(inertia, gradient, u * 0.5, u)


# --- the energy identity, checked against a numerical divergence -----------
def _grid(nt=81, nr=81, moving=True):
    """A front that moves outward while the flow carries an azimuthal component."""
    t = np.linspace(0.0, 4.0, nt)
    r = np.linspace(0.5, 20.0, nr)
    T, R = np.meshgrid(t, r, indexing="ij")
    dt, dr = t[1] - t[0], r[1] - r[0]

    shift = 0.5 * T if moving else 0.0
    beta_r = 0.35 * np.tanh((R - 10.0 - shift) / 2.5)
    beta_phi = 0.25 * np.exp(-(((R - 8.0) / 6.0) ** 2))
    gamma = 1.0 / np.sqrt(1.0 - (beta_r**2 + beta_phi**2))
    zeros = np.zeros_like(T)
    u = np.stack([gamma, gamma * beta_r, zeros, gamma * beta_phi], axis=-1)

    fields = {
        "rho": 1.0 + 0.3 * np.tanh((R - 10.0 - shift) / 2.5),
        "p_t": 0.4 * (1.0 + 0.3 * np.tanh((R - 10.0 - shift) / 2.5)),
        "delta": 0.5 * np.tanh((R - 12.0) / 3.0),
        "pi": 0.1 / np.cosh((R - 11.0) / 3.0) ** 2,
    }
    return t, r, T, R, dt, dr, u, fields


def _kinematics(u, dt, dr):
    """Convective derivative, expansion, shear and the s-stretching, all by finite differences."""
    covariant = u.copy()
    covariant[..., 0] *= -1.0            # u_mu = (-u^0, u^1, u^2, u^3)

    d_u = np.zeros(u.shape[:2] + (4, 4))
    d_u[:, :, 0, :] = np.gradient(covariant, dt, axis=0)
    d_u[:, :, 1, :] = np.gradient(covariant, dr, axis=1)

    theta = np.gradient(u[..., 0], dt, axis=0) + np.gradient(u[..., 1], dr, axis=1)
    return d_u, theta


def test_energy_identity_matches_the_numerical_divergence():
    """``Q`` from the identity equals ``-u.(grad.T_vac)`` to discretization order."""
    t, r, T, R, dt, dr, u, fields = _grid()
    s = radial_direction(u)
    d_u, theta = _kinematics(u, dt, dr)

    stress = vacuum_stress(fields["rho"], fields["p_t"], fields["delta"],
                           u=u, s=s, pi_amplitude=fields["pi"])
    numerical = -_dot(u, flat_divergence(stress, dt, dr))

    h_up = np.diag([-1.0, 1.0, 1.0, 1.0]) + _outer(u, u)
    # The shear lives in the mixed projector h^alpha{}_mu h^beta{}_nu, because d_u acts on a
    # covector; using the all-up h^{mu nu} here would mix the two index positions.
    h_mixed = np.eye(4) + _outer(u, np.einsum("...j,jk->...k", u, np.diag([-1.0, 1.0, 1.0, 1.0])))
    projected = np.einsum("...am,...bn,...ab->...mn", h_mixed, h_mixed, d_u)
    sigma = 0.5 * (projected + np.swapaxes(projected, -1, -2)) - (theta / 3.0)[..., None, None] * h_up
    pi_tensor = fields["pi"][..., None, None] * (_outer(s, s) - h_up / 3.0)

    d_rho = np.zeros(u.shape[:2] + (4,))
    d_rho[:, :, 0] = np.gradient(fields["rho"], dt, axis=0)
    d_rho[:, :, 1] = np.gradient(fields["rho"], dr, axis=1)

    predicted = (
        np.einsum("...a,...a->...", u, d_rho)
        + (fields["rho"] + fields["p_t"]) * theta
        + fields["delta"] * np.einsum("...a,...b,...ab->...", s, s, d_u)
        + np.einsum("...ab,...ab->...", pi_tensor, sigma)
    )

    core = (slice(6, -6), slice(6, -6))
    scale = np.max(np.abs(numerical[core]))
    assert scale > 1e-3                       # the configuration is not trivially zero
    relative = np.max(np.abs((numerical - predicted)[core])) / scale
    assert relative < 2e-2, f"identity violated at {relative:.2e} of the signal"

    # The same sum, computed through the documented interface, agrees term by term.
    for point in ((20, 20), (30, 40), (50, 55)):
        assert np.isclose(
            vacuum_energy_exchange(
                rho_dot=float(np.einsum("...a,...a->...", u, d_rho)[point]),
                rho=float(fields["rho"][point]),
                p_t=float(fields["p_t"][point]),
                theta=float(theta[point]),
                delta=float(fields["delta"][point]),
                direction_stretch=float(
                    np.einsum("...a,...b,...ab->...", s, s, d_u)[point]
                ),
                pi_sigma=float(np.einsum("...ab,...ab->...", pi_tensor, sigma)[point]),
            ),
            float(predicted[point]),
            rtol=1e-9,
            atol=1e-14,
        )
