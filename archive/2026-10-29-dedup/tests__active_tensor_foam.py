"""
Active Tensor Foam Node Model (QVEX Framework)
Canonical implementation of the non-negative foam friction and vacuum density.
"""

import numpy as np

class ActiveTensorFoamNode:
    def __init__(self, rho_0=0.045, r_core=4.5, lambda_coh=12.0):
        self.rho_0 = rho_0          # Vacuum density amplitude
        self.r_core = r_core        # Core radius (kpc)
        self.lambda_coh = lambda_coh  # Coherence length (kpc)

    def vacuum_density(self, r):
        """Lorentzian vacuum density profile."""
        return self.rho_0 / (1.0 + (r / self.r_core) ** 2)

    def screening(self, r):
        """Screening function for active foam."""
        return 1.0 / (1.0 + r / self.lambda_coh)

    def foam_friction(self, r):
        """
        Non-negative nonlinear foam friction profile.
        F_friction(r) = max(0, 1.1 * ln(1 + r/r_c) * (1 - exp(-r/lambda_coh)) - 0.12 * max(0, r - r_c))
        """
        r = np.asarray(r, dtype=float)
        base_activation = 1.1 * np.log1p(r / self.r_core) * (1.0 - np.exp(-r / self.lambda_coh))
        penalty = 0.12 * np.maximum(0.0, r - self.r_core)
        return np.maximum(0.0, base_activation - penalty)

    def compute(self, r_grid):
        """Compute all field quantities on a radial grid."""
        r = np.asarray(r_grid, dtype=float)
        return {
            "r": r,
            "rho": self.vacuum_density(r),
            "screening": self.screening(r),
            "foam_friction": self.foam_friction(r)
        }
