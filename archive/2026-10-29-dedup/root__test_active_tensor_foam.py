"""Unit tests for active_tensor_foam.py using pytest."""

import numpy as np
import pytest
from nodes.active_tensor_foam import ActiveTensorFoamNode


def test_default_parameters():
  node = ActiveTensorFoamNode()
  assert node.lambda_coh == 12.0
  assert node.r_core == 4.5
  assert node.rho_0 == 0.045


def test_vacuum_popcorn_density_at_zero():
  node = ActiveTensorFoamNode()
  assert np.isclose(node.vacuum_density(0.0), node.rho_0)


def test_nonlinear_foam_friction_core():
  node = ActiveTensorFoamNode()
  # At r = 0, friction should be 0
  assert np.isclose(node.foam_friction(0.0), 0.0)


def test_stress_tensor_proxy_positivity():
  node = ActiveTensorFoamNode()
  r_vals = np.linspace(0, 30, 50)
  result = node.compute(r_vals)
  friction = result["foam_friction"]
  # Friction proxy must remain non-negative across physical radii (O1 Guard)
  assert np.all(friction >= 0.0)
