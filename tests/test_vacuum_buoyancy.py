import numpy as np

from nodes.active_tensor_foam import ActiveTensorFoamSolver

def test_buoyancy_positive_and_finite():
    solver = ActiveTensorFoamSolver()
    r_grid = np.linspace(0.1, 30.0, 100)
    buoyancy = solver.vacuum_archimedes_buoyancy(r_grid, g_local=9.81, c_light=1.0)
    
    assert len(buoyancy) == len(r_grid)
    assert np.all(np.isfinite(buoyancy))
    assert np.all(buoyancy >= 0.0)

def test_buoyancy_vanishes_when_homogeneous():
    solver = ActiveTensorFoamSolver()
    r_grid = np.array([1.0, 5.0, 10.0])
    buoyancy_null = solver.vacuum_archimedes_buoyancy(r_grid, background_ratio=1.0)
    
    assert np.allclose(buoyancy_null, 0.0)

def test_core_asymptotics():
    solver = ActiveTensorFoamSolver()
    b_near = solver.vacuum_archimedes_buoyancy(0.5)
    b_far = solver.vacuum_archimedes_buoyancy(50.0)
    
    assert b_near > b_far
