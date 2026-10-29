from pathlib import Path
import sys
import numpy as np

# Allow direct execution from any working directory.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from nodes.active_tensor_foam import ActiveTensorFoamNode

def main():
    print("=" * 60)
    print("Initializing NGC 6503 Active Tensor Foam Fit Pipeline...")
    print("=" * 60)

    node = ActiveTensorFoamNode()
    r_grid = np.linspace(0.1, 30.0, 60)
    results = node.compute(r_grid)

    gamma = results["exchange_coefficient"]
    force = results["exchange_friction"]
    power = node.exchange_power(r_grid, relative_velocity=1.0)
    min_raw_friction = float(np.min(node.foam_friction_raw(r_grid)))
    mean_density = np.mean(results["rho"])

    print(f"Radial grid points evaluated: {len(r_grid)}")
    print(f"Exchange coefficient Gamma (>= 0 by construction): {np.min(gamma):.4f} .. {np.max(gamma):.4f}")
    print(f"Max |exchange force| at du = 1: {np.max(np.abs(force)):.4f}")
    print(f"Min raw Node 02 friction (documented, changes sign): {min_raw_friction:.4e}")
    print(f"Mean vacuum density: {mean_density:.4f}")
    
    assert np.all(gamma >= 0.0), "Violation: the exchange coefficient went negative!"
    assert np.all(force <= 0.0), "Violation: the exchange force no longer opposes relative motion!"
    assert np.all(power <= 0.0), "Violation: the exchange is injecting energy!"
    print("Constraint Check: exchange current is dissipative at every radius [PASS].")
    print("Pipeline execution completed successfully (status: unverified but self-consistent).")
    print("=" * 60)

if __name__ == "__main__":
    main()
