# NGC 6503 SPARC reproducibility bundle

Run with `python reproduce_ngc6503.py` (Python 3.8+, standard library only). Keep `RC_MRT.mrt` beside the script. It checks exactly 31 rows, reports Vobs−Vbar residuals and chi-square, repeats the computation for determinism, and verifies CSV readback.

Important scope: the source QVEX node does not state a mapping from density/friction to circular velocity. Therefore the residuals and reduced chi-square are for the published fixed baryonic Vbar reference, **not** a QVEX velocity prediction or fit. Degrees of freedom are N−0 because no parameters were fitted.

The O1 friction clamp (`max(0, F_raw)`) has been removed. The production drag law is the
Node 05 causal exchange current in `nodes_analysis/05_Vacuum_Exchange_Current_Node.md`,
whose coefficient is non-negative by construction. This script computes neither the foam
fields nor the drag.
