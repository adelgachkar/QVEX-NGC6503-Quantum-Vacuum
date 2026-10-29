"""Extract the canonical equations from ``nodes_analysis`` and evaluate them.

The point of this module is that the Markdown stays the single source of truth for
the model. Tests call :func:`parse_equation` on a node document, evaluate the result
with the implementation's own parameters, and compare it to the implementation. If a
document and the code drift apart, the comparison fails and the LaTeX translation
raises instead of silently skipping.

Only the constructs actually used by the canonical node equations are translated:
``\\frac``, ``^{...}``, ``\\text{...}``, ``\\left``/``\\right``, ``\\exp``, ``\\ln``,
``\\max``, ``\\cdot`` and implicit multiplication. Anything else raises ``ValueError``
rather than evaluating to something wrong.
"""

from pathlib import Path
import re

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "nodes_analysis"

_MACROS = {
    "\\rho_0": "rho_0",
    "\\lambda_{coh}": "lambda_coh",
    "\\gamma_{vac}": "gamma_vac",
    "\\Gamma": "Gamma",
    "\\exp": "exp",
    "\\ln": "log",
    "\\max": "maximum",
    "\\cdot": "*",
    "\\left": "",
    "\\right": "",
}

_NAMESPACE = {
    "exp": np.exp,
    "log": np.log,
    "maximum": np.maximum,
}


def read_node_markdown(filename):
    path = ANALYSIS_DIR / filename
    if not path.is_file():
        raise FileNotFoundError(f"Could not find node file: {path}")
    return path.read_text(encoding="utf-8")


def _matching_brace(text, open_index):
    depth = 0
    for index in range(open_index, len(text)):
        char = text[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index
    raise ValueError(f"unbalanced braces in LaTeX: {text[open_index:]!r}")


def _block_contents(text, opening, closing):
    start = text.index(opening) + len(opening)
    end = text.index(closing, start)
    return text[start:end].strip()


def extract_rhs(content, opening="$$", closing="$$"):
    """Return the right-hand side of the first display equation in ``content``."""
    body = _block_contents(content, opening, closing)
    if "=" not in body:
        raise ValueError(f"no equation found between {opening!r} and {closing!r}")
    return body.split("=", 1)[1].strip()


def latex_to_python(expression):
    """Translate the canonical LaTeX subset into a Python expression."""
    text = expression

    # \text{...} carries no mathematical meaning here.
    while "\\text{" in text:
        index = text.index("\\text{")
        brace = index + len("\\text")
        end = _matching_brace(text, brace)
        text = text[:index] + text[brace + 1:end] + text[end + 1:]

    # \frac{A}{B} -> ((A)/(B))
    while "\\frac" in text:
        index = text.index("\\frac")
        num_open = text.index("{", index)
        num_close = _matching_brace(text, num_open)
        den_open = text.index("{", num_close + 1)
        den_close = _matching_brace(text, den_open)
        numerator = text[num_open + 1:num_close]
        denominator = text[den_open + 1:den_close]
        text = (
            text[:index]
            + f"(({numerator})/({denominator}))"
            + text[den_close + 1:]
        )

    # x^{...} -> x**(...)
    while "^{" in text:
        index = text.index("^{")
        end = _matching_brace(text, index + 1)
        inner = text[index + 2:end]
        text = text[:index] + f"**({inner})" + text[end + 1:]

    for macro, replacement in _MACROS.items():
        text = text.replace(macro, replacement)

    # Remaining subscripts such as R_{cutoff} are plain identifiers in Python.
    while "_{" in text:
        index = text.index("_{")
        end = _matching_brace(text, index + 1)
        text = text[:index] + "_" + text[index + 2:end] + text[end + 1:]

    if "\\" in text:
        raise ValueError(f"untranslated LaTeX command in {text!r}")

    # Implicit multiplication: ")(" and "x(" and "1.1 log(...)".
    text = re.sub(r"\)(?=\s*[\(\w])", ")*", text)
    text = re.sub(r"(?<=[\w\)])\s+(?=[A-Za-z\(])", "*", text)
    text = text.replace(" ", "")
    return text


def parse_equation(filename, opening="$$", closing="$$"):
    """Return a Python source string for the first display equation of a node."""
    return latex_to_python(extract_rhs(read_node_markdown(filename), opening, closing))


def evaluate(expression, solver, r):
    """Evaluate a parsed equation with the solver's documented parameters."""
    namespace = dict(_NAMESPACE)
    namespace.update(
        {
            "r": np.asarray(r, dtype=float),
            "rho_0": solver.rho_0,
            "lambda_coh": solver.lambda_coh,
            "gamma_vac": solver.gamma_vac,
            "R_cutoff": solver.R_cutoff,
            "Gamma_0": solver.friction_scale,
            "r_core": solver.r_core,
        }
    )
    return eval(expression, {"__builtins__": {}}, namespace)  # noqa: S307 - no builtins
