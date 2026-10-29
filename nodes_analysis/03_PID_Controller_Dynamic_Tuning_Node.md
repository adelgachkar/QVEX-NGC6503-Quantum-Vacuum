---
title: "Node 03: PID Controller Dynamic Tuning"
author: "Adel Gachkar"
project: "QVEX-NGC6503-Quantum-Vacuum"
date: "2026-10-26"
tags: [PID_Control, FeedbackLoop, Stability, NGC6503]
---

# Node 03: PID Controller Dynamic Tuning

## Overview

This node implements a Proportional-Integral-Derivative (PID) feedback control loop to dynamically adjust the vacuum coupling parameter $\eta(r)$ along the radial coordinate.

## Control Law and Parameters

The feedback control signal $\eta(r)$ is governed by the tracking error $e(r) = v_{\text{obs}}(r) - v_{\text{model}}(r)$:

$$\eta(r) = K_P e(r) + K_I \int_{0}^{r} e(r') dr' + K_D \frac{de}{dr}$$

> **Notational note:** [[QVEX_NGC6503_Technical_Reference]] §3 writes the same law in the time-domain form $\eta(t)$; the radial form above is canonical for the spatial SPARC application.

### Optimized Gain Matrix:

- Proportional Gain ($K_P$): `0.92`
- Integral Gain ($K_I$): `0.35`
- Derivative Gain ($K_D$): `0.18`

## Stability Analysis

The selected gains guarantee zero steady-state error at the galaxy edge ($r = 11.8 \text{ kpc}$) and maintain phase margin stability across all SPARC radial sampling nodes.
