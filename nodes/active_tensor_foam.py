"""
nodes/active_tensor_foam.py
Active Tensor Foam and Vacuum Buoyancy Model for NGC 6503.

The density and friction laws here are the ones documented in ``nodes_analysis``:

* Node 01 — ``rho_popcorn(r) = rho_0 (1 + r/lambda_coh)^(-gamma_vac) exp(-r/R_cutoff)``
* Node 02 — ``F_friction(r) = 1.1 ln(1 + r/4.5) (1 - exp(-r/12.0)) - 0.12 max(0, r - 4.5)``
* Node 05 — ``Gamma(r) = Gamma_0 ln(1 + r/r_core) (1 - exp(-r/lambda_coh))``

``tests/test_ngc6503_vacuum_nodes.py`` parses those equations out of the Markdown
and compares them numerically against this module, so the document is the source of
truth: changing either side alone makes the suite fail.

The Node 02 expression is kept verbatim in :meth:`foam_friction_raw` because it is what
the document says, and it changes sign near r = 12.11 kpc. It is **not** used as the
production drag law: that role belongs to the Node 05 exchange current, whose coefficient
is non-negative by construction and therefore cannot produce anti-damping at any radius.
The former ``max(0, raw)`` clamp is gone, so nothing has to be hidden.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, Union


def _nonnegative_radius(r: Union[float, np.ndarray]) -> np.ndarray:
    """Validate a radius instead of silently reinterpreting it (Annex A, A.3.1.2)."""
    r_arr = np.asarray(r, dtype=float)
    if np.any(~np.isfinite(r_arr)):
        raise ValueError("radius must be finite")
    if np.any(r_arr < 0.0):
        raise ValueError("radius must be non-negative")
    return r_arr


@dataclass
class ActiveTensorFoamSolver:
    # --- Node 01: vacuum popcorn density -------------------------------------
    rho_0: float = 7.02e-30          # Node 01 amplitude, in g/cm^3
    lambda_coh: float = 12.0         # Node 01 coherence length, in kpc
    gamma_vac: float = 1.0           # Node 01 power-law falloff index
    R_cutoff: float = 50.0           # Node 01 exponential cutoff, in kpc
    density_unit: str = "g/cm^3"

    # --- Node 02: nonlinear foam friction ------------------------------------
    friction_scale: float = 1.1      # Node 02 logarithmic activation amplitude
    friction_penalty: float = 0.12   # Node 02 outer-disk penalty coefficient
    r_core: float = 4.5              # Node 02 core-boundary transition, in kpc

    # --- auxiliary helpers (not part of Nodes 01/02) -------------------------
    alpha: float = 0.8               # screening exponent
    G: float = 4.30091e-6            # kpc (km/s)^2 / Msun

    @property
    def r_s(self) -> float:
        return self.r_core

    def vacuum_density(self, r: Union[float, np.ndarray]) -> np.ndarray:
        """Node 01 canonical profile, in ``density_unit``."""
        r_arr = np.asarray(r, dtype=float)
        return (
            self.rho_0
            * (1.0 + r_arr / self.lambda_coh) ** (-self.gamma_vac)
            * np.exp(-r_arr / self.R_cutoff)
        )

    def screening(self, r: Union[float, np.ndarray]) -> np.ndarray:
        """Auxiliary screening factor: S(0) = 1, S(lambda_coh) = 0.5, S(inf) -> 0."""
        r_arr = np.asarray(r, dtype=float)
        return 1.0 / (1.0 + (r_arr / self.lambda_coh) ** self.alpha)

    def foam_friction_raw(self, r: Union[float, np.ndarray]) -> np.ndarray:
        """Node 02 canonical friction exactly as documented, sign included.

        The documented expression changes sign at r ~ 12.11 kpc and is therefore
        anti-damping over the outer disk. This method keeps that raw behaviour so it
        can be compared with the document and regression-tested; use
        ``foam_friction`` for the non-negative quantity used downstream.
        """
        r_arr = np.asarray(r, dtype=float)
        activation = self.friction_scale * np.log1p(r_arr / self.r_core) * (
            1.0 - np.exp(-r_arr / self.lambda_coh)
        )
        penalty = self.friction_penalty * np.maximum(0.0, r_arr - self.r_core)
        return activation - penalty

    def exchange_coefficient(self, r: Union[float, np.ndarray]) -> np.ndarray:
        """Node 05 exchange coefficient ``Gamma(r)``, non-negative by construction.

        This is the Node 02 activation profile without the linear penalty. Both factors
        are non-negative for ``r >= 0``, so ``Gamma`` cannot change sign and the sign of
        the resulting force is carried by the relative velocity alone. The penalty term
        is what made the documented expression anti-damping beyond 12.11 kpc; it stays
        visible in :meth:`foam_friction_raw` and out of the production law.
        """
        r_arr = _nonnegative_radius(r)
        return self.friction_scale * np.log1p(r_arr / self.r_core) * (
            1.0 - np.exp(-r_arr / self.lambda_coh)
        )

    def exchange_friction(
        self,
        r: Union[float, np.ndarray],
        relative_velocity: Union[float, np.ndarray] = 1.0,
    ) -> np.ndarray:
        """``F = -Gamma * dv``, the non-relativistic limit of ``J^mu = Gamma (u_v^mu - u_m^mu)``.

        The force opposes the *relative* motion at every radius: it depends on a velocity
        difference rather than on an absolute velocity or a velocity gradient, so it is
        Galilean invariant and cannot turn into a drive.
        """
        velocity = np.asarray(relative_velocity, dtype=float)
        return -self.exchange_coefficient(r) * velocity

    def exchange_power(
        self,
        r: Union[float, np.ndarray],
        relative_velocity: Union[float, np.ndarray] = 1.0,
    ) -> np.ndarray:
        """``u_m . J = -Gamma * dv^2 <= 0``: the exchange can only remove energy."""
        velocity = np.asarray(relative_velocity, dtype=float)
        return -self.exchange_coefficient(r) * velocity ** 2

    def exchange_balance(
        self,
        r: Union[float, np.ndarray],
        relative_velocity: Union[float, np.ndarray] = 1.0,
    ) -> Dict[str, np.ndarray]:
        """Matter source and vacuum sink of one conserved pair (Annex A, A.3.3.2)."""
        matter_source = self.exchange_friction(r, relative_velocity)
        return {
            'matter_source': matter_source,
            'vacuum_sink': -matter_source,
        }

    def vacuum_archimedes_buoyancy(
        self,
        r: Union[float, np.ndarray],
        g_local: float = 9.81,
        c_light: float = 1.0,
        background_ratio: float = 0.0,
    ) -> np.ndarray:
        r_arr = np.asarray(r, dtype=float)
        rho = self.vacuum_density(r_arr)
        eps = 1e-5
        rho_plus = self.vacuum_density(r_arr + eps)
        rho_minus = self.vacuum_density(np.maximum(0.0, r_arr - eps))
        drho_dr = (rho_plus - rho_minus) / (2.0 * eps)
        buoyancy = - (1.0 - background_ratio) * (rho / self.rho_0) * drho_dr * g_local
        return np.asarray(buoyancy, dtype=float)

    def effective_enclosed_mass(self, r: Union[float, np.ndarray]) -> np.ndarray:
        """Auxiliary rotation-curve helper — not part of Node 01 or Node 02.

        Known limitation: it feeds ``rho_0`` in ``density_unit`` into a solver whose
        ``G`` is in kpc/Msun units, i.e. the density is not converted to
        Msun/kpc^3 (7.02e-30 g/cm^3 ~ 1.04e2 Msun/kpc^3). No test binds this helper.
        """
        r_arr = np.asarray(r, dtype=float)
        return 4.0 * np.pi * self.rho_0 * (self.r_core ** 3) * (
            np.log(1.0 + r_arr / self.r_core) - (r_arr / (self.r_core + r_arr))
        )

    def circular_velocity(self, r: Union[float, np.ndarray]) -> np.ndarray:
        r_arr = np.asarray(r, dtype=float)
        m_enc = self.effective_enclosed_mass(r_arr)
        safe_r = np.where(r_arr == 0, 1e-9, r_arr)
        v_sq = self.G * m_enc / safe_r
        v_sq = np.clip(v_sq, 0.0, None)
        return np.sqrt(v_sq) * self.screening(r_arr)

    def compute(self, r: Union[list, np.ndarray]) -> Dict[str, np.ndarray]:
        r_arr = np.asarray(r, dtype=float)
        return {
            'r': r_arr,
            'rho': self.vacuum_density(r_arr),
            'screening': self.screening(r_arr),
            'exchange_coefficient': self.exchange_coefficient(r_arr),
            'exchange_friction': self.exchange_friction(r_arr),
        }


class ActiveTensorFoamNode(ActiveTensorFoamSolver):
    pass
