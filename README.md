"""
===============================================================================
Ocenographic data processing script
File: oceanographic_data_processing.py
===============================================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Purpose
-------
This appendix documents the Python-based computational workflow developed to
derive near-bed hydrodynamic parameters, including bed shear stresses,
wave–current interaction metrics, turbulent mixing-layer thickness (δ_mix),
and statistical descriptors for each sampling interval (C0–C7). The analysis
is based on synchronized near-bed current and wave records.

Scientific Basis
----------------
The workflow follows classical formulations for sediment transport and
wave–current interaction, including those of Bagnold (1956),
Soulsby and Smallman (1986), and Soulsby and Clarke (2005). These formulations
are used to compute current-induced shear stress (τ_c), wave-induced shear
stress (τ_w), mean shear stress (τ_m), and the maximum combined bed shear
stress (τ_max).

The turbulent mixing-layer thickness δ_mix is computed using the
parameterization of Harris and Wiberg (1997):

    δ_mix = 0.07 · (τ_max,int − τ_cr) + 6 · D50

where:
    τ_max,int = maximum τ_max observed within each sampling interval (C0–C7),
    τ_cr      = critical shear stress for incipient motion,
    D50       = median grain size of the native sediment.

Two distinct critical shear-stress thresholds (τ_cr) are employed:
    1. A native-sediment τ_cr (based on D50), used in δ_mix computations and in
       evaluating the potential for sediment mobilization under combined flow.
    2. A tracer-based τ_cr (based on the finer tracer D50), used exclusively for
       interpretive assessment of whether hydrodynamic conditions were capable
       of mobilizing tracer particles. This value does not influence δ_mix.

The framework ensures full computational reproducibility of all hydrodynamic
forcing metrics presented in Sections 4 and 5 of the main text.

Input Data Assumptions
---------------------
The input Excel file must contain synchronized, burst-averaged current and
wave measurements with the following mandatory fields:
    ve, vn   : near-bed current velocity components (m s⁻¹)
    D50      : median native sediment grain size (m)
    z        : current measurement elevation above bed (m)
    h        : water depth (m)
    Tz, Hs   : zero-crossing wave period (s) and significant wave height (m)
    dir_w    : wave direction (degrees)
    Time     : timestamp (datetime format)
    tau_cr   : critical shear stress for native sediment (N m⁻²)

Records lacking any mandatory fields are excluded from processing.

Implementation Details
----------------------
- Language:   Python 3.10+
- Libraries: pandas, numpy, openpyxl, argparse, logging
- Input:     Excel workbook containing burst-averaged current and wave parameters
- Output:    Structured Excel workbook with computed stresses, δ_mix values,
             and interval-based summary statistics
- Structure: Modular functions implementing physical formulations, data
             processing, statistical aggregation, and export routines
"""
