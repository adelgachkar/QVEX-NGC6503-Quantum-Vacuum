"""One vacuum density, one unit.

Node 01 quotes the amplitude in ``g/cm^3``, Annex A2 bounds it in ``kg/m^3`` and the
rotation-curve helper needs ``Msun/kpc^3``. Those are three representations of one number,
a factor of 1e3 and 1e31 apart. Before this module existed, the number was copied between
domains as a bare float, so a value in the wrong unit looked exactly like a value in the
right one.

These tests parse the unit table out of *both* documents, convert, and compare. Editing
one column, or restating the amplitude in a fourth unit, turns the suite red with the two
readings printed side by side.
"""

from pathlib import Path
import re

import numpy as np
import pytest

from nodes.active_tensor_foam import (
    ActiveTensorFoamSolver,
    CANONICAL_DENSITY_UNIT,
    DENSITY_UNITS,
    convert_density,
)

ROOT = Path(__file__).resolve().parents[1]
NODE_01 = ROOT / "nodes_analysis" / "01_Vacuum_Popcorn_Density_Node.md"
ANNEX_A2 = ROOT / "nodes_analysis" / "Annex_A2_Vacuum_Popcorn_Boundary_Conditions.md"

#: Documented amplitudes are rounded to three significant figures in the Msun column.
MSUN_RTOL = 6e-3


def documented_table(path):
    """Return ``{symbol: {unit: value}}`` parsed from the canonical unit table.

    The table has one row per parameter and one column per unit; the column heading names
    the unit, so the parsed values carry their unit with them instead of being bare floats.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = None
    for index, line in enumerate(lines):
        if line.strip().startswith("|") and "kg/m^3" in line and "g/cm^3" in line:
            header_index = index
            break
    if header_index is None:
        raise AssertionError(f"{path.name} has no canonical unit table")

    units = [cell.strip() for cell in lines[header_index].strip().strip("|").split("|")]
    unit_columns = {}
    for column, heading in enumerate(units):
        match = re.search(r"(kg/m\^3|g/cm\^3|Msun/kpc\^3)", heading)
        if match:
            unit_columns[match.group(1)] = column

    table = {}
    for line in lines[header_index + 1:]:
        if not line.strip().startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if set("".join(cells)) <= set("-: "):
            continue
        symbol = cells[1]
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", symbol):
            continue
        table[symbol] = {
            unit: float(cells[column]) for unit, column in unit_columns.items()
        }
    return table


@pytest.fixture(scope="module")
def node_table():
    return documented_table(NODE_01)


@pytest.fixture(scope="module")
def annex_table():
    return documented_table(ANNEX_A2)


def test_every_document_publishes_the_same_amplitude_in_every_unit(node_table, annex_table):
    """The red test: two documents, one value.

    Both tables are read column by column and every column is converted into the canonical
    unit before comparison. A number written in the wrong column cannot pass, because the
    comparison is between quantities, not between strings.
    """
    assert "rho_0" in node_table, "Node 01 lost its canonical unit table"
    assert "rho_0" in annex_table, "Annex A2 lost its canonical unit table"

    readings = {}
    for name, table in (("Node 01", node_table), ("Annex A2", annex_table)):
        for unit, value in table["rho_0"].items():
            readings.setdefault(unit, {})[name] = value

    mismatches = []
    for unit, by_document in readings.items():
        values = list(by_document.values())
        if len(values) < 2:
            mismatches.append(f"{unit}: only one document states this column {by_document}")
            continue
        canonical = [float(convert_density(value, unit, CANONICAL_DENSITY_UNIT)) for value in values]
        rtol = MSUN_RTOL if unit == "Msun/kpc^3" else 1e-12
        # atol=0 is not decoration: the amplitudes are ~1e-27, far below numpy's default
        # atol of 1e-8, so a default comparison would call a factor-1000 error "close".
        if not np.allclose(canonical, canonical[0], rtol=rtol, atol=0.0):
            mismatches.append(
                f"{unit}: {by_document} -> {canonical} {CANONICAL_DENSITY_UNIT}"
            )
    assert not mismatches, (
        "the two documents disagree about the vacuum density once converted to "
        f"{CANONICAL_DENSITY_UNIT}: {mismatches}"
    )


def test_each_table_is_internally_consistent(node_table, annex_table):
    """Within one document, every unit column must describe the same density."""
    for name, table in (("Node 01", node_table), ("Annex A2", annex_table)):
        row = table["rho_0"]
        canonical = {
            unit: float(convert_density(value, unit, CANONICAL_DENSITY_UNIT))
            for unit, value in row.items()
        }
        reference = canonical["kg/m^3"]
        for unit, value in canonical.items():
            rtol = MSUN_RTOL if unit == "Msun/kpc^3" else 1e-12
            assert np.isclose(value, reference, rtol=rtol, atol=0.0), (
                f"{name}: rho_0 in {unit} converts to {value} {CANONICAL_DENSITY_UNIT}, "
                f"but the kg/m^3 column says {reference}"
            )


def test_the_code_default_is_the_documented_amplitude(node_table):
    solver = ActiveTensorFoamSolver()
    documented_kg = node_table["rho_0"]["kg/m^3"]
    assert solver.density_unit in DENSITY_UNITS
    assert np.isclose(solver.rho_0_kg_m3, documented_kg, rtol=1e-12, atol=0.0)
    assert np.isclose(
        float(convert_density(solver.rho_0, solver.density_unit, "g/cm^3")),
        node_table["rho_0"]["g/cm^3"],
        rtol=1e-12,
        atol=0.0,
    )
    assert np.isclose(
        solver.rho_0_msun_kpc3, node_table["rho_0"]["Msun/kpc^3"], rtol=MSUN_RTOL, atol=0.0
    )


def test_documented_bound_is_satisfied_in_the_canonical_unit(annex_table):
    """A bound in kg/m^3 must be checked against the amplitude in kg/m^3.

    The bound line is located by text, not by a LaTeX regex, so the check survives
    cosmetic edits to the formula while still failing if the unit of the bound changes.
    """
    annex = ANNEX_A2.read_text(encoding="utf-8")
    bound_lines = [
        line for line in annex.splitlines()
        if "rho_0" in line and "10^{" in line and "kg/m" in line and "^3" in line
    ]
    assert bound_lines, "Annex A2 no longer states the rho_0 bound in the canonical unit"

    exponent = int(bound_lines[0].split("10^{", 1)[1].split("}", 1)[0])
    upper = 10.0 ** exponent
    amplitude = annex_table["rho_0"]["kg/m^3"]
    assert 0.0 <= amplitude <= upper, (
        f"rho_0 = {amplitude} kg/m^3 violates the documented bound {upper} kg/m^3"
    )
    # The same value must be a thousand times smaller when expressed in g/cm^3; if the
    # amplitude were ever written in the wrong column, one of these two readings would fail.
    assert np.isclose(
        float(convert_density(amplitude, "kg/m^3", "g/cm^3")),
        annex_table["rho_0"]["g/cm^3"],
        rtol=1e-12,
        atol=0.0,
    )


def test_the_withdrawn_symbols_do_not_come_back():
    """One parameter, one name: r_s / r_c / alpha must not be reintroduced next to them."""
    annex = ANNEX_A2.read_text(encoding="utf-8")

    # The profile is the display equation that carries the exponential cutoff.
    profile = [
        block for block in re.findall(r"\$\$(.+?)\$\$", annex, re.S)
        if r"\rho_{\text{popcorn}}(r) =" in block and r"\exp" in block
    ]
    assert profile, "Annex A2 no longer states the density profile"
    equation = profile[0]

    assert r"\frac{r}{\lambda_{\text{coh}}}" in equation, (
        "Annex A2 must state the profile with Node 01's lambda_coh"
    )
    assert r"R_{\text{cutoff}}" in equation
    assert r"\gamma_{\text{vac}}" in equation
    for withdrawn in (r"{r_s}", r"{r_c}", r"^{-\alpha}"):
        assert withdrawn not in equation, (
            f"Annex A2 reintroduced the withdrawn symbol {withdrawn} in its profile equation"
        )


def test_conversion_is_exact_and_reversible():
    for unit in DENSITY_UNITS:
        value = 3.7e-27
        round_trip = convert_density(
            convert_density(value, "kg/m^3", unit), unit, "kg/m^3"
        )
        assert np.isclose(round_trip, value, rtol=1e-15, atol=0.0)

    # The two units that were confused in this bundle differ by exactly 1000.
    assert np.isclose(convert_density(1.0, "g/cm^3", "kg/m^3"), 1000.0, rtol=1e-15)
    assert np.isclose(
        convert_density(7.02e-30, "g/cm^3", "kg/m^3"), 7.02e-27, rtol=1e-15, atol=0.0
    )


def test_unknown_units_raise_instead_of_being_guessed():
    with pytest.raises(ValueError):
        convert_density(1.0, "g/cm3", "kg/m^3")
    with pytest.raises(ValueError):
        ActiveTensorFoamSolver().vacuum_density(1.0, unit="solar masses")


def test_the_profile_can_only_be_reinterpreted_by_conversion():
    solver = ActiveTensorFoamSolver()
    grid = np.linspace(0.0, 30.0, 25)
    native = solver.vacuum_density(grid)
    converted = solver.vacuum_density(grid, unit="Msun/kpc^3")

    assert np.isclose(native[0], solver.rho_0, rtol=1e-15, atol=0.0)
    assert np.allclose(
        converted,
        convert_density(native, solver.density_unit, "Msun/kpc^3"),
        rtol=1e-14,
        atol=0.0,
    )
    # 1.037e2, not 7.02e-30: the converted profile is not the native array.
    assert not np.allclose(converted, native, rtol=1e-14, atol=0.0)


def test_enclosed_mass_converts_the_density_before_using_msun_gravity():
    """The solver's G is in kpc/Msun, so the density must arrive in Msun/kpc^3."""
    solver = ActiveTensorFoamSolver()
    radii = np.array([1.0, 5.0, 20.0])
    mass = solver.effective_enclosed_mass(radii)

    rho_msun = solver.rho_0_msun_kpc3
    expected = 4.0 * np.pi * rho_msun * solver.r_core ** 3 * (
        np.log1p(radii / solver.r_core) - radii / (solver.r_core + radii)
    )
    assert np.allclose(mass, expected, rtol=1e-14)

    # Independent check: with the documented 1.037e2 Msun/kpc^3 the halo mass inside
    # 1 kpc is 2239 Msun. Feeding the raw g/cm^3 number instead gives 1.5e-28, which is the
    # factor-1.5e31 unit error this conversion exists to prevent.
    assert np.isclose(mass[0], 2239.2, rtol=2e-3), mass[0]
