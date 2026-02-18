"""
====================================================
Sediment mobility parameters and coordinate rotation
File: sediment_mobility_and_rotation.py
====================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Purpose
-------
This workflow computes key sediment mobility parameters following
Soulsby (1997) and applies a shoreline-oriented coordinate rotation
to spatial datasets. All computations are performed using a single
Excel input file containing sediment properties and geographic
coordinates for each sampling point.

Input
-----
Input Excel file:
    INPUT_SEDIMENT_PARAMETERS.xlsx

Each worksheet (e.g., Tracer, C2, C3, C4, C5, C6, …) must contain:
    - x, y      : original spatial coordinates (m)
    - D50       : median sediment grain size (m)
    - rho_sed   : sediment density (kg/m³)
    - rho_sw    : seawater density (kg/m³)
    - nu        : kinematic viscosity (m²/s)
    - g         : gravitational acceleration (m/s²)

Optional columns:
    - Sample    : sample identifier or code
    - Sed_Camp  : campaign name (created automatically if missing)
    - delta_mix : mixing-layer thickness (if provided)

Output
------
OUTPUT_SEDIMENT_MOBILITY_PARAMETERS.xlsx containing:
    * Sediment_Parameters   : D*, θ_cr, τ_cr, u*_cr, w_s, and related metrics
    * Rotated_Coordinates   : original (x, y) and rotated (as, cs) coordinates
    * Campaign_Statistics   : per-sample campaign table (one row per sample)
    * Summary_Statistics    : grouped mean/std (by Local & Sed_Camp)
    * Global_Statistics     : total sample count and global summary statistics

Additional output:
    * COORDINATE_ROTATION_DIAGRAM.jpg (quality-assurance plot)

Notes
-----
Along-shore (as) and cross-shore (cs) coordinates are obtained by
rotating the original (x, y) positions about the tracer origin
(the tracer injection point). The default rotation angle is −22°,
corresponding to a clockwise rotation that aligns the coordinate
system with the local shoreline orientation.
"""
