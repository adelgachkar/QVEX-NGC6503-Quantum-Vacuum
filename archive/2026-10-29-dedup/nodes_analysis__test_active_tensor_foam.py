"""Unit tests for active_tensor_foam.py using pytest."""

import numpy as np
from active_tensor_foam import ActiveTensorFoamSolver, FoamParameters
import pytest


def test_default_parameters():
  solver = ActiveTensorFoamSolver()
  assert solver.params.lambda_coh == 12.0
  assert solver.params.r_core == 4.5


def test_vacuum_popcorn_density_at_zero():
  solver = ActiveTensorFoamSolver()
  rho_0 = solver.params.rho_0
  assert np.isclose(solver.vacuum_popcorn_density(0.0), rho_0)


def test_nonlinear_foam_friction_core():
  solver = ActiveTensorFoamSolver()
  # At r = 0, friction should be 0
  assert np.isclose(solver.nonlinear_foam_friction(0.0), 0.0)


def test_stress_tensor_proxy_positivity():
  solver = ActiveTensorFoamSolver()
  r_vals = np.linspace(0, 30, 50)
  pi_eff = solver.compute_stress_tensor_proxy(r_vals)
  # Stress proxy must remain non-negative across physical radii
  assert np.all(pi_eff >= 0.0)
