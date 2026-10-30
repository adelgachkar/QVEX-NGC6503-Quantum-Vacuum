"""Invariants of the non-equilibrium front stress (Node 06).

Each test states a property that must hold for *any* parameters, not a golden number, so
it cannot be satisfied by tuning a constant. The two that carry the physics argument are
``test_leading_and_trailing_fronts_carry_opposite_stress`` (the equilibrium ``Gamma(r)``
law cannot produce this) and ``test_memory_leaves_the_front_out_of_equilibrium`` (the
front is non-equilibrium exactly when the Deborah number is order one or larger).
"""

import numpy as np
import pytest

from nodes.non_equilibrium_stress import (
    FrontStressState,
    anisotropy_pressure,
    causal_speed,
    crossing_time,
    deborah_number,
    entropy_production,
    equilibrium_stress,
    exchange_decomposition,
    four_velocity,
    front_shear,
    front_summary,
    geometry,
    is_orthogonal_to_flow,
    projector,
    radial_shear,
    relaxed_stress,
    shell_velocity,
    spatial_direction,
    stress_tensor,
    trace,
    validate_transport,
)

SHELL = dict(inner_radius=10.0, outer_radius=20.0, width=0.5, velocity=0.1)


# --- tensor structure -------------------------------------------------------
@pytest.mark.parametrize("delta", [0.0, 1.0, -2.5, 7.25])
def test_stress_is_orthogonal_to_flow_traceless_and_symmetric(delta):
    """u.Pi = 0, Pi^mu_mu = 0, Pi^mu nu = Pi^nu mu — for every anisotropy."""
    tensor = stress_tensor(delta)

    assert is_orthogonal_to_flow(tensor)
    assert np.isclose(trace(tensor), 0.0, atol=1e-12)
    assert np.allclose(tensor, tensor.T)
    # Longitudinal vs transverse split of a traceless radial stress.
    assert np.isclose(tensor[1, 1], 2.0 * delta / 3.0)
    assert np.isclose(tensor[2, 2], -delta / 3.0)
    assert np.isclose(tensor[3, 3], -delta / 3.0)


def test_projector_has_the_expected_rank_and_signature():
    """h = diag(0,1,1,1): it removes the time direction and keeps three spatial ones."""
    h = projector()
    assert np.allclose(h, np.diag([0.0, 1.0, 1.0, 1.0]))
    assert np.isclose(trace(h), 3.0)
    u = four_velocity()
    assert np.allclose(h @ u, 0.0)


def test_the_geometry_rejects_a_direction_that_is_not_orthogonal_to_the_flow():
    """A radial direction with a time component is not a spacelike rest-frame axis."""
    # Unit spacelike, but tilted into the time direction: normalized yet not orthogonal.
    tilted = spatial_direction() + 0.5 * four_velocity()
    tilted = tilted / np.sqrt(float(tilted @ geometry() @ tilted))
    assert np.isclose(float(tilted @ geometry() @ tilted), 1.0)
    with pytest.raises(ValueError, match="orthogonal to the flow"):
        stress_tensor(1.0, s=tilted)

    with pytest.raises(ValueError, match="timelike"):
        stress_tensor(1.0, u=np.array([1.0, 0.1, 0.0, 0.0]))


def test_anisotropy_round_trips_through_the_tensor():
    for delta in (-3.0, 0.5, 4.0):
        recovered = anisotropy_pressure(np.asarray(stress_tensor(delta))[1, 1])
        assert np.isclose(float(recovered), delta)


# --- transport --------------------------------------------------------------
def test_equilibrium_stress_opposes_the_shear_and_is_antisymmetric_in_it():
    positive = equilibrium_stress(eta=2.0, shear=0.3)
    negative = equilibrium_stress(eta=2.0, shear=-0.3)

    assert np.isclose(float(positive), -2.0 * 2.0 * 0.3)
    assert np.isclose(float(positive), -float(negative))
    assert np.sign(positive) != np.sign(0.3)


def test_relaxation_follows_the_analytic_exponential():
    tau, pi_eq = 3.0, -1.4
    assert np.isclose(float(relaxed_stress(pi_eq, 0.0, tau)), 0.0)
    assert np.isclose(float(relaxed_stress(pi_eq, tau, tau)), pi_eq * (1 - np.exp(-1)))
    assert np.isclose(float(relaxed_stress(pi_eq, 12 * tau, tau)), pi_eq, rtol=1e-4)
    # Memory: the sign of the initial state matters, which the equilibrium law cannot see.
    assert float(relaxed_stress(pi_eq, 0.0, tau, initial_value=5.0)) == 5.0


def test_entropy_production_is_non_negative_for_the_dissipative_branch():
    """-(1/2) Pi : sigma = eta sigma^2 >= 0, at both fronts and for either sign of sigma."""
    for sigma in (-0.4, 0.4):
        pi = float(equilibrium_stress(eta=1.5, shear=sigma))
        assert float(entropy_production(pi, sigma)) > 0.0
        # Pi : sigma itself is negative; reporting that as entropy production would
        # reverse the arrow of time for a dissipative front.
        assert pi * sigma < 0.0


# --- front kinematics -------------------------------------------------------
def test_one_shell_has_two_fronts_with_opposite_shear():
    """The whole point: a single velocity field gives +sigma and -sigma."""
    shear = front_shear(**SHELL)
    assert shear["leading"] < 0.0 < shear["trailing"]
    assert np.isclose(shear["leading"], -shear["trailing"])

    # The analytic peak is what the field actually carries at each front's centre.
    assert np.isclose(radial_shear(SHELL["outer_radius"], **SHELL), shear["leading"])
    assert np.isclose(radial_shear(SHELL["inner_radius"], **SHELL), shear["trailing"])


def test_shell_velocity_is_zero_outside_and_flat_inside():
    assert np.isclose(shell_velocity(0.0, **SHELL), 0.0, atol=1e-9)
    assert np.isclose(shell_velocity(50.0, **SHELL), 0.0, atol=1e-9)
    middle = 0.5 * (SHELL["inner_radius"] + SHELL["outer_radius"])
    assert np.isclose(shell_velocity(middle, **SHELL), SHELL["velocity"], rtol=1e-6)

    with pytest.raises(ValueError):
        shell_velocity(middle, 20.0, 10.0, 0.5, 0.1)


def test_crossing_time_and_deborah_number():
    assert np.isclose(crossing_time(0.5, 0.1), 5.0)
    assert np.isclose(deborah_number(2.0, 0.5, 0.1), 0.4)
    with pytest.raises(ValueError):
        crossing_time(0.0, 1.0)


def test_memory_leaves_the_front_out_of_equilibrium():
    """De << 1 -> the equilibrium law is right; De >> 1 -> the stress has not arrived."""
    fast = front_summary(eta=0.3, tau=0.2, inertia=10.0, **SHELL)
    slow = front_summary(eta=0.3, tau=20.0, inertia=10.0, **SHELL)

    assert fast["deborah_number"] < 0.1 and not fast["is_non_equilibrium"]
    assert slow["deborah_number"] > 3.0 and slow["is_non_equilibrium"]

    fast_ratio = fast["fronts"]["leading"]["relaxation_ratio"]
    slow_ratio = slow["fronts"]["leading"]["relaxation_ratio"]
    assert fast_ratio > 0.99
    assert slow_ratio < 0.3
    assert slow_ratio > 0.0

    # The limits are Monotone in De: more memory, less stress actually delivered.
    ordered = [
        front_summary(eta=0.3, tau=tau, inertia=10.0, **SHELL)["fronts"]["leading"][
            "relaxation_ratio"
        ]
        for tau in (0.1, 1.0, 5.0, 50.0)
    ]
    assert ordered == sorted(ordered, reverse=True)


def test_both_fronts_are_dissipative_and_opposite_in_sign():
    summary = front_summary(eta=0.3, tau=20.0, inertia=10.0, **SHELL)
    leading = summary["fronts"]["leading"]
    trailing = summary["fronts"]["trailing"]

    assert leading["entropy_production"] > 0.0
    assert trailing["entropy_production"] > 0.0
    assert leading["shear"] * trailing["shear"] < 0.0
    assert leading["stress_at_crossing"] * trailing["stress_at_crossing"] < 0.0
    assert summary["stress_opposes_shear"] is True


# --- causality --------------------------------------------------------------
def test_causality_is_enforced_rather_than_assumed():
    """tau -> 0 is the parabolic limit: the stress would propagate instantly."""
    assert np.isclose(causal_speed(eta=9.0, tau=1.0, inertia=10.0), np.sqrt(0.9))
    with pytest.raises(ValueError, match="causality"):
        validate_transport(eta=20.0, tau=1.0, inertia=10.0)  # v_c = sqrt(2) > 1
    with pytest.raises(ValueError, match="tau_Pi"):
        validate_transport(eta=1.0, tau=0.0, inertia=10.0)
    with pytest.raises(ValueError, match="anti-damping"):
        validate_transport(eta=-1.0, tau=1.0, inertia=10.0)
    with pytest.raises(ValueError):
        front_summary(eta=1e3, tau=1.0, inertia=1.0, **SHELL)


# --- interaction vector -----------------------------------------------------
def test_exchange_vector_splits_into_energy_and_force_and_conserves_the_pair():
    split = exchange_decomposition(Q_energy=2.5, force_radial=-0.75)
    vector = split["vector"]
    u, g = four_velocity(), geometry()

    assert np.isclose(float((g @ u) @ vector), -2.5)          # u.Q = -Q
    assert np.isclose(vector[1], -0.75)                       # f^alpha = h^alpha_nu Q^nu
    assert np.allclose(split["matter_source"] + split["vacuum_sink"], 0.0)
    assert np.isclose(vector[0], 2.5)                         # the energy component is Q u^0


def test_front_state_reports_its_own_invariants():
    state = FrontStressState(shear=-0.0667, stress=0.4)
    tensor = state.tensor

    assert is_orthogonal_to_flow(tensor)
    assert np.isclose(trace(tensor), 0.0, atol=1e-12)
    assert np.isclose(state.anisotropy, 0.6)
    assert state.is_dissipative
