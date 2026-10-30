"""What the ``phi-r`` shear channel can and cannot do.

Three claims are tested separately, because they are usually run together:

1. ``Pi = Delta chi`` (Nodes 06 and 07) has **no** ``r-phi`` component, so a longitudinal
   anisotropy cannot torque anything. Reaching ``Pi^{r phi} != 0`` requires the *shear*
   part ``-2 eta sigma``, which needs differential rotation.
2. With that shear the transport is real: angular momentum flows outward for a falling
   rotation curve, and the torque density is ``eta v / r^2`` for a flat one.
3. The transport is **not** centripetal support and is not free: the torque must be paid
   for at the rate ``eta v^2 / r^2``, and the same stress dissipates
   ``eta v^2 / (2 r^2)``.

The shear formula is derived here from the coordinate-basis covariant derivative rather
than asserted: ``test_orbital_shear_matches_the_coordinate_basis_derivative`` recomputes it
with finite differences and the connection term ``Gamma^phi_{r phi} = 1/r``.
"""

import numpy as np
import pytest

from nodes.angular_momentum import (
    LIGHT_SPEED_KPC_PER_MYR,
    angular_momentum_flux,
    decay_rate_needed,
    dissipation,
    flat_curve_power,
    flat_curve_stress,
    flat_curve_torque_density,
    orbital_shear,
    power_budget,
    torque_density,
    torque_power,
    viscous_stress,
)
from nodes.non_equilibrium_stress import stress_tensor, trace


# --- the shear --------------------------------------------------------------
def test_orbital_shear_matches_the_coordinate_basis_derivative():
    """Recompute sigma_{r phi} from d_r(r^2 Omega) - (2/r) r^2 Omega and compare."""
    radius = np.linspace(2.0, 20.0, 2000)
    omega = 0.1175 / radius**1.5           # Keplerian, so the shear is large
    u_phi = radius**2 * omega              # the coordinate component, gamma = 1

    # sigma_{r phi} = 1/2[d_r u_phi - (2/r) u_phi]; the orthonormal component divides by r.
    d_u_phi = np.gradient(u_phi, radius)
    sigma_coordinate = 0.5 * (d_u_phi - 2.0 * u_phi / radius)
    sigma_orthonormal = sigma_coordinate / radius

    analytic = orbital_shear(radius[1:-1], omega[1:-1], np.gradient(omega, radius)[1:-1])
    # Second-order finite differences and a second-order gradient: agreement to ~1e-4.
    assert np.allclose(sigma_orthonormal[1:-1], analytic, rtol=1e-3, atol=1e-9)


def test_solid_body_rotation_carries_no_shear():
    """The channel needs differential rotation; a rigid rotor has none."""
    assert np.isclose(orbital_shear(10.0, 0.01, 0.0), 0.0)


def test_flat_and_keplerian_shears():
    v_flat, radius = 0.1175, 10.0
    omega_flat = v_flat / radius
    sigma_flat = orbital_shear(radius, omega_flat, -v_flat / radius**2)
    assert np.isclose(sigma_flat, -v_flat / (2.0 * radius))

    omega_kep = np.sqrt(1.0) / radius**1.5
    sigma_kep = orbital_shear(radius, omega_kep, -1.5 * omega_kep / radius)
    assert np.isclose(sigma_kep, -0.75 * omega_kep)


# --- the two candidate channels --------------------------------------------
@pytest.mark.parametrize("delta", [-3.0, 0.0, 2.5])
def test_longitudinal_anisotropy_has_no_r_phi_component(delta):
    """``Pi = Delta chi`` cannot torque: chi^{r phi} = 0 for a radial s, identically."""
    tensor = stress_tensor(delta)
    # Exactly zero, not merely small: no parameter choice makes this component non-zero.
    assert tensor[1, 3] == 0.0 and tensor[3, 1] == 0.0
    assert np.isclose(trace(tensor), 0.0, atol=1e-12)


def test_the_shear_part_is_what_carries_the_torque():
    eta, shear = 195.0, -5.875e-3
    stress = viscous_stress(eta, shear)
    assert stress > 0.0
    assert np.isclose(stress, -2.0 * eta * shear)

    # Falling rotation curve -> positive stress -> outward transport.
    assert angular_momentum_flux(10.0, stress) > 0.0
    # And reversing the sign of dOmega/dr reverses the transport direction.
    assert viscous_stress(eta, -shear) < 0.0


def test_negative_viscosity_is_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        viscous_stress(-1.0, 0.1)
    with pytest.raises(ValueError, match="positive"):
        orbital_shear(0.0, 0.1, -0.01)


# --- the balance ------------------------------------------------------------
def test_torque_density_is_the_divergence_of_the_angular_momentum_flux():
    """(1/r^2) d_r(r^2 T) computed from the flux, versus the closed form."""
    radius = np.linspace(3.0, 20.0, 300)
    eta, v_flat = 195.0, 0.1175
    stress = eta * v_flat / radius
    flux = radius**2 * stress
    numerical = np.gradient(flux, radius) / radius**2

    closed = np.array([torque_density(r, s, -eta * v_flat / r**2)
                       for r, s in zip(radius, stress)])
    assert np.allclose(numerical[2:-2], closed[2:-2], rtol=1e-3, atol=1e-9)
    # Same numbers as the direct flat-curve expression.
    assert np.allclose(closed, [flat_curve_torque_density(eta, v_flat, r) for r in radius])


def test_flat_curve_needs_a_torque_that_scales_with_eta_and_with_v_squared():
    """A rigid flat curve is not torque-free: the shear is -v/(2r) everywhere."""
    radius, v_flat = 10.0, 0.1175
    assert np.isclose(flat_curve_stress(100.0, v_flat, radius),
                      100.0 * v_flat / radius)
    # Halving the viscosity halves the torque; doubling v doubles the stress but quadruples
    # the power.
    torque = flat_curve_torque_density(195.0, v_flat, radius)
    assert np.isclose(flat_curve_torque_density(97.5, v_flat, radius), torque / 2.0)
    assert np.isclose(flat_curve_power(195.0, 2 * v_flat, radius),
                      4.0 * flat_curve_power(195.0, v_flat, radius))


# --- the energy price -------------------------------------------------------
def test_dissipation_is_half_the_torque_work_for_a_flat_curve():
    """The known disk relation, checked rather than quoted: P_torque = 2 * D."""
    radius, v_flat, eta = 10.0, 0.1175, 195.0
    sigma = -v_flat / (2.0 * radius)

    assert np.isclose(dissipation(eta, sigma), eta * v_flat**2 / (2.0 * radius**2))
    assert np.isclose(flat_curve_power(eta, v_flat, radius), 2.0 * dissipation(eta, sigma))
    # ... and the torque's own work rate equals flat_curve_power.
    torque = flat_curve_torque_density(eta, v_flat, radius)
    assert np.isclose(torque_power(radius, torque, v_flat / radius),
                      flat_curve_power(eta, v_flat, radius))


def test_power_budget_reports_the_deficit_rather_than_assuming_support():
    radius, v_flat, eta = 10.0, 0.1175, 195.0
    sigma = -v_flat / (2.0 * radius)
    torque = flat_curve_torque_density(eta, v_flat, radius)
    needed = flat_curve_power(eta, v_flat, radius) + dissipation(eta, sigma)

    unsupported = power_budget(0.0, torque, radius, v_flat / radius, eta, sigma)
    assert unsupported["sustained"] is False
    assert np.isclose(unsupported["required_power"], needed)
    assert np.isclose(unsupported["net"], -needed)

    exact = power_budget(needed, torque, radius, v_flat / radius, eta, sigma)
    assert exact["sustained"] is True and abs(exact["net"]) < 1e-12

    generous = power_budget(2.0 * needed, torque, radius, v_flat / radius, eta, sigma)
    assert generous["sustained"] is True and generous["net"] > 0.0


def test_the_price_of_the_torque_is_a_transient_vacuum():
    """Quantitative verdict: express the required power in orbital periods, with c explicit.

    No threshold is assumed. The assertion records what the numbers say for the canonical
    vacuum amplitude and the largest shear the causality bound of Node 06 permits: the
    vacuum would be spent in *less than one rotation*, so it cannot be a static halo
    sustaining the curve.
    """
    radius, v_flat, eta = 10.0, 0.1175, 195.0
    sigma = -v_flat / (2.0 * radius)
    torque = flat_curve_torque_density(eta, v_flat, radius)
    needed = flat_curve_power(eta, v_flat, radius) + dissipation(eta, sigma)

    rho_v = 103.7                                   # Msun/kpc^3, from Node 01
    rest_energy = rho_v * LIGHT_SPEED_KPC_PER_MYR**2
    budget = decay_rate_needed(needed, rest_energy)

    orbital_period = 2.0 * np.pi * radius / v_flat  # Myr
    assert budget["efolding_myr"] < orbital_period
    assert np.isclose(budget["efolding_myr"], 1.0 / budget["rate_per_myr"])

    # The rate scales with eta, so a smaller eta is cheaper -- but the torque is smaller by
    # the same factor, so the mechanism stops mattering rather than becoming free.
    cheap = decay_rate_needed(needed / 10.0, rest_energy)
    assert np.isclose(cheap["rate_per_myr"] * 10.0, budget["rate_per_myr"])
    # No power, no decay needed.
    assert decay_rate_needed(0.0, rest_energy)["efolding_myr"] == float("inf")
