"""
================================================================
Sediment textural analysis — Folk & Ward (1957) graphical method
File: sediment_textural_analysis_folk_ward.py
================================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Description
-----------
This script computes grain-size statistical parameters following the
Folk and Ward (1957) graphical (percentile-based) method using cumulative
grain-size distributions expressed on the Krumbein φ scale.

Calculated parameters include:
    - Mean grain size (Mz)
    - Sorting (σI)
    - Skewness (SkI)
    - Kurtosis (KG)

Textural classifications (grain size, sorting, skewness, and kurtosis)
are assigned using standard Folk and Ward threshold criteria.

Input
-----
Excel file containing:
    • Column 1: sample name or identifier
    • Columns 2 to n: grain-size classes (φ or mm) with associated
      frequency values (%)

Units are automatically detected unless explicitly specified via the
command-line argument --units {auto | phi | mm}. Grain-size classes are
sorted internally from coarse to fine prior to analysis.

Output
------
An Excel workbook containing:
    • Textural_Parameters_phi — percentiles, Mz, σI, SkI, KG, and associated
      textural classes
    • Percentiles_mm          — percentile values converted to millimetres
    • QA                      — quality-assurance diagnostics including
                                 monotonicity checks, percentile validity,
                                 zero-spread flags, and summary statistics
    • Phi_Scale_Reference     — reference table for the Krumbein φ scale
    • Classification_Criteria — Folk and Ward classification ranges and
                                 governing formulas

This script provides a reproducible and transparent workflow for the
textural analysis of native sediments and tracer materials using the
Folk and Ward (1957) graphical method.
"""
