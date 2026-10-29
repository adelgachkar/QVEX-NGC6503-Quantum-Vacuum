import numpy as np

from nodes.active_tensor_foam import ActiveTensorFoamNode


def test_active_tensor_foam_initialization():
    node = ActiveTensorFoamNode()
    assert (node.r_core, node.lambda_coh) == (4.5, 12.0)
    assert node.density_unit == "g/cm^3"


def test_density_and_screening():
    node = ActiveTensorFoamNode()
    assert np.isclose(node.vacuum_density(0.0), node.rho_0)
    assert np.isclose(node.screening(0.0), 1.0)
    assert np.isclose(node.screening(12.0), 0.5)


def test_exchange_coefficient_is_nonnegative_and_zero_at_origin():
    values = ActiveTensorFoamNode().exchange_coefficient(np.linspace(0, 50, 100))
    assert np.all(values >= 0.0)
    assert values[0] == 0.0


def test_raw_friction_is_exposed_for_document_comparison():
    node = ActiveTensorFoamNode()
    grid = np.linspace(0, 50, 100)
    raw = node.foam_friction_raw(grid)

    assert raw[0] == 0.0
    assert np.any(raw < 0.0), "the documented expression must keep its outer-disk sign change"
    assert not hasattr(node, "foam_friction"), "the clamp was replaced by the Node 05 exchange law"


def test_compute_returns_documented_keys_and_shapes():
    result = ActiveTensorFoamNode().compute([0.0, 5.0, 10.0])
    assert set(result) == {
        'r', 'rho', 'screening', 'exchange_coefficient', 'exchange_friction'
    }
    assert all(len(value) == 3 for value in result.values())
