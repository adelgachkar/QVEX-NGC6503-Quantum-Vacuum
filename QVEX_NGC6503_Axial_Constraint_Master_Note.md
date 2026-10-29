---
title: "NGC 6503 Rotation Analysis Using an Axial Constraint in the QVEX Framework"
author: "Adel Gachkar"
project: "QVEX-NGC6503-Quantum-Vacuum"
document_type: "Independent Research Note"
date: "2026-10-26"
tags:
  - QVEX
  - NGC6503
  - QuantumVacuum
  - PIDControl
  - ObsidianMaster
---

## NGC 6503 Rotation Analysis Using an Axial Constraint in the QVEX Framework

**Author:** Adel Gachkar  
**Status:** Independent research note, maintained separately from the main code repository.

---

### 1. Overview and Related Notes

This document serves as the central note for the QVEX analysis of the galaxy NGC 6503. Detailed calculations and model-specific discussions are maintained in the following related notes:

- [[01_Vacuum_Popcorn_Density_Node|Vacuum fluctuation density]]
- [[02_Nonlinear_Foam_Friction_Node|Nonlinear foam friction]]
- [[03_PID_Controller_Dynamic_Tuning_Node|Dynamic tuning of the PID controller]]
- [[04_Dipole_Axis_Constraint_Node|Dipole-axis constraint and central jets]]

Related documents:

- [[NGC6503_QVEX_Complete_Analysis|Complete analysis report]]
- [[QVEX_NGC6503_Technical_Reference|Technical reference note]]
- [[QVEX_NGC6503_Review_Report_v1|Technical review report (v1)]]

---

### 2. Summary of Model Parameters

Below is the complete and unified summary of the model parameters, friction formulation, and visual reference paths:

```text
=== PID Controller Gains ===
KP = 0.92
KI = 0.35
KD = 0.18

=== Nonlinear Foam-Friction Function ===
F_friction(r) = 1.1 * ln(1.0 + r / 4.5) * (1.0 - exp(-r / 12.0)) - 0.12 * max(0, r - 4.5)

=== Reference Plot ===
outputs/plots/sparc_ngc6503_fit.png
```
