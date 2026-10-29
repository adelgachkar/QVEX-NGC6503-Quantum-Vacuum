# Node 02: Nonlinear Quantum Foam Friction

This node models dissipative, resistive drag on test masses moving through spacetime quantum foam at galactic peripheries (specifically verified for NGC 6503).

## 1. Mathematical Formulation & Implementation

The canonical empirically calibrated radial friction force, matched to Technical Reference §2 and the `test_nonlinear_foam_friction_node_o1_regression` test, is formulated as:

\[
F_{\text{friction}}(r) = 1.1 \ln\left(1 + \frac{r}{4.5}\right) \left(1 - \exp\left(-\frac{r}{12.0}\right)\right) - 0.12 \max(0, r - 4.5)
\]

### Key Physical Constraints & Calibration:

- **Core-Boundary Transition:** Occurs at $r = 4.5\,\text{kpc}$, marking the onset of outer-disk stabilization.
- **Asymptotic Damping:** The exponential damping term $\left(1 - \exp(-r/12.0)\right)$ ensures that unphysical outer-edge velocity divergences are suppressed, guaranteeing stable asymptotic convergence toward the characteristic plateau velocity ($v_\infty \approx 120\,\text{km/s}$).
- **Non-Negativity / O1 Status (corrected 2026-10-29):** The linear penalty term $-0.12 \max(0, r-4.5)$ does **not** keep this expression non-negative. It makes $F_{\text{friction}}$ change sign near $r \approx 12.11\,\text{kpc}$, so the expression is anti-damping beyond that radius; the earlier statement that the penalty ensures physical non-negativity was wrong. This radial expression is kept here and in `foam_friction_raw()` as the documented law for regression comparison. The production drag law is the [[05_Vacuum_Exchange_Current_Node|Node 05 causal exchange current]], whose coefficient is non-negative by construction and needs no clamp.

## 2. Conceptual Velocity-Gradient Precursor

As a theoretical precursor, the underlying microscopic damping is expressed via the velocity-gradient form:

\[
F_{\text{foam}}(v, r) = -\gamma_{\text{foam}} v(r) \left|\frac{\partial v}{\partial r}\right|^\alpha, \qquad \alpha = 1.5
\]

Here, $\gamma_{\text{foam}}$ is calibrated against the standard MOND acceleration scale $a_0 = 1.15 \times 10^{-10}\,\text{m/s}^2$. The radial formulation above serves as the rigorously tested, production-ready effective model in the QVEX pipeline.
