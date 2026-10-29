---
title: "Technical Reference & Analytical Note: Independent Validation of NGC 6503 via QVEX Vacuum-Foam Framework"
author: "Adel Gachkar, Islamic Azad University, Urmia"
orcid: "0009-0006-7713-6004"
email: "adelgachkar@gmail.com"
reference_id: "QVEX-REF-NGC6503-2026-V1"
repository: "github.com/adelgachkar/QVEX_v0.1.0"
---

## 1. Introduction and Epistemic Separation

This technical reference document provides an independent, citable summary of the application of the Quantum Vacuum Exchange (QVEX) engine to the SPARC database rotation curve of NGC 6503 (Lelli, McGaugh & Schombert 2016). All SPARC fitting routines, non-linear vacuum friction derivations, and PID stabilization nodes are documented herein as a standalone engineering report to maintain strict version integrity of the core release (v0.1.0).

## 2. Theoretical Foundations: Vacuum Weighting and Buoyancy Formalism

The QVEX framework interprets galactic vacuum foam as a dense medium interacting with a gravitational field. To provide a rigorous theoretical anchor for this interaction, we adopt the weighting methodology established by Shevchenko & Shevrin (2016) regarding Archimedean forces in Casimir systems.

### 2.1 Covariant Conservation and Vacuum Energy Weighting

The core interaction of vacuum energy (foam density) with the metric is defined by the covariant conservation of the total energy-momentum tensor:
$$\nabla_\mu T^{\mu\nu} = 0$$

Following Shevchenko & Shevrin, the effective gravitational weight ($f$) of the vacuum-foam region in a weak gravitational field is expressed via integrated pressure and energy density weighting:
$$f = \frac{1}{2} \int d^3x \sqrt{-g_{00}(z)} \, \frac{\partial g_{00}(z)}{\partial z} \left[ \langle T^{00} \rangle_2 - \langle T^{00} \rangle_1 \right]$$

This confirms that the vacuum energy density ($\langle T^{00} \rangle$) directly contributes to the effective gravitational force via quantum vacuum buoyancy. This justifies our treatment of vacuum-foam density ($\rho_{\text{foam}}$) as a buoyant medium that responds to galactic gravitational gradients, eliminating the necessity for classical dark matter halos.

### 2.2 Epistemic Bridge to Galactic Metrics

The effective dipole potential ($V_{\text{eff}}$) used in our NGC 6503 analysis is now epistemically harmonized with this weighting procedure. The jet-axis constraint (derived from the foam-friction node) is a local manifestation of the vacuum buoyancy effect in the presence of an asymmetric field gradient, as necessitated by the equivalence principle for quantum-field media.

Moreover, rotation-velocity flattening is coupled to a scale-dependent vacuum coherence length $\lambda_{\text{coh}} = 12.0\text{ kpc}$. The non-linear vacuum friction function is defined as:
$$F_{\text{friction}}(r) = 1.1 \cdot \ln\left(1.0 + \frac{r}{4.5}\right) \cdot \left(1.0 - \exp\left(-\frac{r}{12.0}\right)\right) - 0.12 \cdot \max\left(0, r - 4.5\right)$$

## 3. PID Feedback Controller & Boundary Stabilization

A closed-loop PID controller regulates the quantum vacuum density gradient $\eta(r)$ to eliminate discretization errors and radial oscillation artifacts at outer boundaries ($r > 15\text{ kpc}$):
$$\eta(t) = K_P \cdot e(t) + K_I \cdot \int e(\tau)d\tau + K_D \cdot \frac{de}{dt}$$

- **$K_P = 0.92$**:

Instantaneous restoring stiffness against gravitational potential gradients.

- **$K_I = 0.35$**: Accumulating inward vacuum flux to eliminate steady-state velocity drift.
- **$K_D = 0.18$**: Viscous quantum foam damping along the rotational axis.

## 4. Central Jet Dipole Constraint & Causal Gauge Limits

Disk stability and axial confinement are maintained via opposing relativistic jets establishing a potential well along the $\hat{z}$ axis:
$$V_{\text{eff}}(r,\theta) = -\frac{GM(r)}{r} + \alpha_{\text{dipole}} \cdot P_{\text{jet}} \cdot \frac{\cos^2(\theta)}{r^3}$$

Standard paradigms (NFW halos and MOND) encounter tensions with localized boundary oscillations and phase-coherent vacuum backreaction. QVEX resolves these via electrodynamic vacuum exchange driven by baryonic currents.

## 5. Verification and Availability

Numerical outputs have been fully verified against SPARC datasets, with corresponding visualization plots stored in `outputs/plots/sparc_ngc6503_fit.png`. Node derivations are cross-referenced in:

- `01_Vacuum_Popcorn_Density_Node`
- `02_Nonlinear_Foam_Friction_Node`
- `03_PID_Controller_Dynamic_Tuning_Node`
- `04_Dipole_Axis_Constraint_Node`
