"""Canonical tests binding the nodes_analysis documents to the executable model.

Node 01 and Node 02 are the source of truth: their display equations are parsed out
of the Markdown by :mod:`tests.node_parser` and evaluated with the implementation's
own parameters. If a document and ``nodes.active_tensor_foam`` disagree — in form, in
a parameter value, or in an inserted literal — these tests fail.

Nodes 03 and 04 have no executable counterpart yet, so they are only checked for the
content they are supposed to carry.
"""

from pathlib import Path
import re

import numpy as np
import pytest

from nodes.active_tensor_foam import ActiveTensorFoamNode, ActiveTensorFoamSolver
from tests.node_parser import evaluate, parse_equation, read_node_markdown

NODE_01 = "01_Vacuum_Popcorn_Density_Node.md"
NODE_02 = "02_Nonlinear_Foam_Friction_Node.md"
NODE_03 = "03_PID_Controller_Dynamic_Tuning_Node.md"
NODE_04 = "04_Dipole_Axis_Constraint_Node.md"
NODE_05 = "05_Vacuum_Exchange_Current_Node.md"
ANNEX_A = "Annex_A_Vacuum_Popcorn_Physical_Constraints.md"
ANNEX_A2 = "Annex_A2_Vacuum_Popcorn_Boundary_Conditions.md"

NODE_FILES = [NODE_01, NODE_02, NODE_03, NODE_04, NODE_05, ANNEX_A, ANNEX_A2]

# rtol used for every document-vs-code comparison. The two sides are the same
# expression evaluated twice, so agreement should be exact to float rounding.
RTOL = 1e-12

RADII = np.linspace(0.0, 50.0, 201)


def parse_documented_scalars(content):
    """Parse ``name = <number>`` parameter lines out of a node document."""
    found = {}
    for key, value in re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([0-9.]+)", content, re.M):
        found.setdefault(key, float(value))
    return found


def parse_documented_scientific(content, key):
    """Parse ``key = <mantissa> x 10^<exponent>`` (times sign or ASCII x)."""
    match = re.search(
        rf"^\s*{key}\s*=\s*([0-9.]+)\s*[×xX*]\s*10\^\(\s*(-?[0-9]+)\s*\)", content, re.M
    ) or re.search(rf"^\s*{key}\s*=\s*([0-9.]+)\s*[×xX*]\s*10\^(-?[0-9]+)", content, re.M)
    if not match:
        raise AssertionError(f"{key} is not documented in scientific notation")
    return float(match.group(1)) * 10.0 ** int(match.group(2))


def test_node_markdown_documents_exist_and_are_nonempty():
    for filename in NODE_FILES:
        assert read_node_markdown(filename).strip(), f"{filename} is empty"


# --- Node 01: vacuum popcorn density ---------------------------------------------

def test_node_01_equation_parses_without_untranslated_latex():
    expression = parse_equation(NODE_01)
    assert "rho_0" in expression
    assert "lambda_coh" in expression
    assert "gamma_vac" in expression
    assert "R_cutoff" in expression


def test_node_01_documented_parameters_match_the_implementation():
    documented = parse_documented_scalars(read_node_markdown(NODE_01))
    node = ActiveTensorFoamSolver()

    assert documented["lambda_coh"] == node.lambda_coh
    assert documented["gamma_vac"] == node.gamma_vac
    assert documented["R_cutoff"] == node.R_cutoff
    assert np.isclose(
        parse_documented_scientific(read_node_markdown(NODE_01), "rho_0"),
        node.rho_0,
        rtol=1e-12,
    ), "documented rho_0 and the implementation amplitude disagree"


def test_node_01_density_matches_the_documented_equation():
    """Divergence test: the Markdown equation must reproduce the code exactly."""
    node = ActiveTensorFoamNode()
    documented = evaluate(parse_equation(NODE_01), node, RADII)
    implemented = node.vacuum_density(RADII)

    assert documented.shape == implemented.shape
    assert np.allclose(documented, implemented, rtol=RTOL, atol=0.0), (
        "Node 01 equation and vacuum_density() disagree: "
        f"max relative deviation {np.max(np.abs(documented - implemented) / np.abs(documented)):.3e}"
    )


def test_node_01_density_profile_behaves_as_documented():
    node = ActiveTensorFoamNode()
    rho = node.vacuum_density(RADII)

    assert np.all(np.isfinite(rho))
    assert np.all(rho >= 0.0)
    assert np.isclose(rho[0], node.rho_0)
    assert np.all(np.diff(rho) <= 0.0), "density must not grow with radius"


# --- Node 02: nonlinear foam friction -------------------------------------------

def test_node_02_equation_parses_without_untranslated_latex():
    expression = parse_equation(NODE_02, r"\[", r"\]")
    for literal in ("1.1", "4.5", "12.0", "0.12"):
        assert literal in expression, f"Node 02 equation lost the constant {literal}"


def test_node_02_raw_friction_matches_the_documented_equation():
    """Divergence test: the Markdown expression must reproduce the raw code exactly."""
    node = ActiveTensorFoamNode()
    documented = evaluate(parse_equation(NODE_02, r"\[", r"\]"), node, RADII)
    implemented = node.foam_friction_raw(RADII)

    assert np.allclose(documented, implemented, rtol=RTOL, atol=1e-15), (
        "Node 02 equation and foam_friction_raw() disagree: "
        f"max absolute deviation {np.max(np.abs(documented - implemented)):.3e}"
    )


def test_node_02_clamp_is_removed_and_the_production_law_never_reverses():
    """The clamp must not come back, and the law that replaced it holds one sign."""
    node = ActiveTensorFoamNode()

    assert not hasattr(node, "foam_friction"), (
        "the max(0, F_raw) clamp must stay removed: positivity is a property of the "
        "Node 05 exchange coefficient, not a post-processing step"
    )

    force = node.exchange_friction(RADII, relative_velocity=1.0)
    assert np.all(force <= 0.0), "the exchange force must oppose relative motion everywhere"
    assert np.all(force[1:] < 0.0), "the exchange force must be strictly damping beyond r = 0"


def test_node_02_sign_change_is_located_and_documented():
    """The outer-disk anti-damping onset must stay at the audited radius."""
    node = ActiveTensorFoamNode()
    grid = np.linspace(0.0, 50.0, 50001)
    raw = node.foam_friction_raw(grid)
    crossing = grid[np.argmax(raw < 0.0)]

    assert 12.0 < crossing < 12.2, f"friction sign change moved to r = {crossing:.2f} kpc"
    assert np.isclose(node.foam_friction_raw(0.0), 0.0)


# --- Nodes 03 and 04: documented-only (no executable counterpart yet) -----------

def test_node_03_documents_the_pid_gains():
    content = read_node_markdown(NODE_03)
    for gain in ("0.92", "0.35", "0.18"):
        assert gain in content, f"Node 03 lost the documented gain {gain}"


def test_node_04_documents_the_dipole_potential():
    content = read_node_markdown(NODE_04)
    assert "V_{\\text{eff}}" in content
    assert "\\alpha_{\\text{dipole}}" in content
    assert "P_{\\text{jet}}" in content


# --- Node 05: causal exchange current -------------------------------------------

def test_node_05_equation_parses_without_untranslated_latex():
    expression = parse_equation(NODE_05)
    assert expression.startswith("Gamma_0*"), expression
    assert "r_core" in expression and "lambda_coh" in expression


def test_node_05_exchange_coefficient_matches_the_documented_equation():
    """Divergence test: the Node 05 equation must reproduce exchange_coefficient()."""
    node = ActiveTensorFoamNode()
    documented = evaluate(parse_equation(NODE_05), node, RADII)
    implemented = node.exchange_coefficient(RADII)

    assert np.allclose(documented, implemented, rtol=RTOL, atol=0.0), (
        "Node 05 equation and exchange_coefficient() disagree: "
        f"max absolute deviation {np.max(np.abs(documented - implemented)):.3e}"
    )


def test_node_05_documented_parameters_match_the_implementation():
    documented = parse_documented_scalars(read_node_markdown(NODE_05))
    node = ActiveTensorFoamNode()

    assert documented["Gamma_0"] == node.friction_scale
    assert documented["r_core"] == node.r_core
    assert documented["lambda_coh"] == node.lambda_coh


def test_node_05_coefficient_is_nonnegative_everywhere():
    node = ActiveTensorFoamNode()
    gamma = node.exchange_coefficient(RADII)

    assert np.all(gamma >= 0.0)
    assert gamma[0] == 0.0
    assert np.all(np.diff(gamma) > 0.0), "the coefficient must grow monotonically outward"


def test_node_05_force_sign_is_set_by_relative_motion_only():
    """This is the structural replacement for the O1 clamp."""
    node = ActiveTensorFoamNode()

    for velocity in (-2.0, -1.0, 1.0, 2.0):
        force = node.exchange_friction(RADII, relative_velocity=velocity)
        assert np.all(np.sign(force[1:]) == -np.sign(velocity)), (
            f"an r-dependent sign reversal reappeared for relative velocity {velocity}"
        )

    # The old law still reverses, which is exactly why it is no longer the production law.
    assert np.any(node.foam_friction_raw(RADII) < 0.0)


def test_node_05_dissipation_never_becomes_a_drive():
    node = ActiveTensorFoamNode()

    for velocity in (-3.0, -1.0, 1.0, 3.0):
        assert np.all(node.exchange_power(RADII, velocity) <= 0.0)


def test_node_05_matter_and_vacuum_terms_sum_to_zero():
    node = ActiveTensorFoamNode()
    balance = node.exchange_balance(RADII, relative_velocity=1.0)

    assert set(balance) == {"matter_source", "vacuum_sink"}
    assert np.all(balance["matter_source"] + balance["vacuum_sink"] == 0.0), (
        "the exchange must conserve the pair exactly, not approximately"
    )


def test_node_05_rejects_a_negative_radius():
    node = ActiveTensorFoamNode()

    with pytest.raises(ValueError):
        node.exchange_coefficient(-1.0)
    with pytest.raises(ValueError):
        node.exchange_friction(np.array([0.0, -0.5]), 1.0)


# --- structural guard ----------------------------------------------------------

def test_only_one_test_module_per_basename_outside_the_archive():
    """Guard the regression this refactor removed: duplicate test basenames."""
    repo_root = Path(__file__).resolve().parents[1]
    basenames = {}
    for path in repo_root.rglob("test_*.py"):
        rel = path.relative_to(repo_root)
        if rel.parts[0] in {".venv", "archive", "__pycache__"}:
            continue
        basenames.setdefault(path.name, []).append(str(rel))
    duplicates = {name: paths for name, paths in basenames.items() if len(paths) > 1}
    assert not duplicates, f"duplicate test module basenames: {duplicates}"
