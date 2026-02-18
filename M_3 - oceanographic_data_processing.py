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


# =============================================================================
# 1. CONFIGURATION
# =============================================================================

import pandas as pd
import numpy as np
import logging
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Dict

# -- Physical constants -------------------------------------------------------
@dataclass
class OceanographicConstants:
    GRAVITY: float = 9.81
    RHO_SEAWATER: float = 1025.0
    ROUGHNESS_FACTOR: float = 2.5
    VON_KARMAN: float = 0.41
    WAVE_PERIOD_FACTOR: float = 1.281  # Tp / Tz ratio

# -- Interval configuration ---------------------------------------------------
@dataclass
class TimeInterval:
    start_date: str
    end_date: str
    sheet_name: str

# -- Processing configuration -------------------------------------------------
@dataclass
class ProcessingConfig:
    input_file: str = 'INPUT_CURRENTS_WAVES_DATA_BASE.xlsx'
    output_file: str = 'OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx'
    tau_max_thresholds: Tuple[float, ...] = (0.000, 0.259, 0.303)

# Global definitions
CONSTANTS = OceanographicConstants()
CONFIG = ProcessingConfig()

TIME_INTERVALS = [
    TimeInterval('2020-08-18 17:50:08', '2020-08-30 23:50:08', 'C0_C2'),
    TimeInterval('2020-08-31 00:05:08', '2020-09-28 23:50:08', 'C2_C3'),
    TimeInterval('2020-09-29 00:05:08', '2020-11-14 23:50:08', 'C3_C4'),
    TimeInterval('2020-11-15 00:05:08', '2021-01-14 23:50:08', 'C4_C5'),
    TimeInterval('2021-01-15 00:05:08', '2021-04-22 23:50:08', 'C5_C6'),
    TimeInterval('2021-04-23 00:05:08', '2021-10-27 23:50:08', 'C6_C7')
]

# Note: we use 'tau_cr' (all lower-case) inside the code.
# If the Excel uses 'Tau_cr' or 'TAU_CR', we normalize later in main().
REQUIRED_COLUMNS = ['ve', 'vn', 'D50', 'z', 'h', 'Tz', 'Hs', 'dir_w', 'Time', 'tau_cr']


# =============================================================================
# 2. PHYSICAL AND OCEANOGRAPHIC FORMULATIONS
# =============================================================================

def calculate_current_magnitude(ve, vn):
    return np.hypot(ve, vn)

def calculate_current_azimuth(ve, vn):
    return (np.degrees(np.arctan2(ve, vn)) + 360) % 360

def calculate_wave_current_angle(cur_dir, wave_dir):
    return np.degrees(np.arctan2(np.sin(np.radians(wave_dir - cur_dir)),
                                 np.cos(np.radians(wave_dir - cur_dir))))

def calculate_nikuradse_roughness(D50):
    return CONSTANTS.ROUGHNESS_FACTOR * D50

def calculate_bed_roughness_length(kn):
    return kn / 30.0

def calculate_natural_scaling_period(h):
    return np.sqrt(h / CONSTANTS.GRAVITY)

def calculate_peak_wave_period(Tz):
    return CONSTANTS.WAVE_PERIOD_FACTOR * Tz

def calculate_period_ratio(Tn, Tz):
    return Tn / Tz

def calculate_wave_correction_factor(t):
    return (6500 + (0.56 + 15.54 * t) ** 6) ** (1 / 6)

def calculate_orbital_velocity_amplitude(Hs, Tn, Wcf, t):
    return (0.25 * Hs) / (Tn * (1 + Wcf * t ** 2) ** 3)

def calculate_max_orbital_velocity(urms):
    return np.sqrt(2) * urms

def calculate_semi_orbital_excursion(uw, Tp):
    return (uw * Tp) / (2 * np.pi)

def calculate_relative_roughness(A, kn):
    return A / kn

def calculate_current_friction_factor(z, z0):
    # Soulsby/Prandtl-von Kármán log-law based bulk friction factor
    return 2 * (CONSTANTS.VON_KARMAN / np.log(z / z0)) ** 2

def calculate_shear_velocity(fc, uz):
    # u_star = u_z * sqrt(fc/2)
    return uz * np.sqrt(fc / 2)

def calculate_current_shear_stress(fc, uz):
    return 0.5 * CONSTANTS.RHO_SEAWATER * fc * uz ** 2

def calculate_wave_friction_factor(r):
    return 0.237 * r ** (-0.52)

def calculate_wave_shear_stress(fwr, uw):
    return 0.5 * CONSTANTS.RHO_SEAWATER * fwr * uw ** 2

def calculate_mean_bed_shear_stress(tau_c, tau_w):
    ratio = tau_w / (tau_c + tau_w)
    return tau_c * (1 + 1.2 * ratio ** 3.2)

def calculate_max_bed_shear_stress(tau_m, tau_w, phi_deg):
    phi = np.radians(phi_deg)
    return np.hypot(tau_m + tau_w * np.cos(phi), tau_w * np.sin(phi))


# =============================================================================
# 3. DATA PROCESSING WORKFLOW
# =============================================================================

def process_oceanographic_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df['Time'] = pd.to_datetime(df['Time'], errors='coerce')

    # Currents
    df['u_z'] = calculate_current_magnitude(df['ve'], df['vn'])
    df['dir_c'] = calculate_current_azimuth(df['ve'], df['vn'])

    # Roughness and current friction factor
    df['kn'] = df['D50'].apply(calculate_nikuradse_roughness)
    df['z0'] = df['kn'].apply(calculate_bed_roughness_length)
    df['fc'] = df.apply(lambda r: calculate_current_friction_factor(r['z'], r['z0']), axis=1)

    # Shear velocity (u_star)
    df['u_star'] = df.apply(lambda r: calculate_shear_velocity(r['fc'], r['u_z']), axis=1)

    # Current shear stress
    df['tau_c'] = df.apply(lambda r: calculate_current_shear_stress(r['fc'], r['u_z']), axis=1)

    # Waves
    df['Tn'] = df['h'].apply(calculate_natural_scaling_period)
    df['Tp'] = df['Tz'].apply(calculate_peak_wave_period)
    df['t']  = df.apply(lambda r: calculate_period_ratio(r['Tn'], r['Tz']), axis=1)
    df['W_cf'] = df['t'].apply(calculate_wave_correction_factor)
    df['urms'] = df.apply(lambda r: calculate_orbital_velocity_amplitude(r['Hs'], r['Tn'], r['W_cf'], r['t']), axis=1)
    df['uw']   = df['urms'].apply(calculate_max_orbital_velocity)
    df['A']    = df.apply(lambda r: calculate_semi_orbital_excursion(r['uw'], r['Tp']), axis=1)
    df['r']    = df.apply(lambda r: calculate_relative_roughness(r['A'], r['kn']), axis=1)
    df['fwr']  = df['r'].apply(calculate_wave_friction_factor)
    df['tau_w'] = df.apply(lambda r: calculate_wave_shear_stress(r['fwr'], r['uw']), axis=1)

    # Coupling and maxima
    df['angle_phi'] = df.apply(lambda r: calculate_wave_current_angle(r['dir_c'], r['dir_w']), axis=1)
    df['tau_m'] = df.apply(lambda r: calculate_mean_bed_shear_stress(r['tau_c'], r['tau_w']), axis=1)
    df['tau_max'] = df.apply(lambda r: calculate_max_bed_shear_stress(r['tau_m'], r['tau_w'], r['angle_phi']), axis=1)

    # -----------------------------------------------------------------
    # Harris & Wiberg (1997) mixing layer thickness δ_mix
    # Using the maximum τ_max in each sampling period (TIME_INTERVALS)
    # δ_mix = 0.07 (τ_max_interval_max − τ_cr) + 6 D50
    # -----------------------------------------------------------------
    # First, create a column to hold the interval-wise τ_max maximum
    df['tau_max_interval'] = np.nan

    for interval in TIME_INTERVALS:
        mask = (
            (df['Time'] >= pd.to_datetime(interval.start_date)) &
            (df['Time'] <= pd.to_datetime(interval.end_date))
        )
        if not mask.any():
            continue

        tau_max_max = df.loc[mask, 'tau_max'].max()
        df.loc[mask, 'tau_max_interval'] = tau_max_max

    # Now compute δ_mix using τ_max_interval (constant within each period)
    # tau_cr is read from Excel (already normalized to lowercase in main()).
    df['delta_mix'] = 0.07 * np.maximum(df['tau_max_interval'] - df['tau_cr'], 0.0) + 6.0 * df['D50']

    return df


# =============================================================================
# 4. STATISTICAL ANALYSIS BY INTERVAL
# =============================================================================

def calculate_interval_statistics(df: pd.DataFrame,
                                  intervals: List[TimeInterval],
                                  thresholds: Tuple[float, ...]) -> Dict:
    """Compute τ_max statistics per time interval and threshold."""
    stats = {'interval_stats': [], 'threshold_stats': []}
    for interval in intervals:
        mask = (
            (df['Time'] >= pd.to_datetime(interval.start_date)) &
            (df['Time'] <= pd.to_datetime(interval.end_date))
        )
        subset = df.loc[mask]
        if subset.empty:
            continue

        row = {
            'Interval': interval.sheet_name,
            'Mean_tau_max': subset['tau_max'].mean(),
            'Max_tau_max': subset['tau_max'].max(),
            'Min_tau_max': subset['tau_max'].min(),
            'Mean_u_star': subset['u_star'].mean(),
            'Max_u_star': subset['u_star'].max(),
            'Min_u_star': subset['u_star'].min(),
            'Count': len(subset)
        }
        stats['interval_stats'].append(row)

        for thr in thresholds:
            sel = subset[subset['tau_max'] > thr]
            if not sel.empty:
                stats['threshold_stats'].append({
                    'Interval': f"{interval.sheet_name}_tau_max>{thr:.3f}",
                    'Threshold': thr,
                    'Count': len(sel),
                    'Max_tau_max': sel['tau_max'].max(),
                    'Mean_u_star_over_thr': sel['u_star'].mean()
                })

    return {k: pd.DataFrame(v) for k, v in stats.items() if v}


# =============================================================================
# 5. OUTPUT AND REPRODUCIBILITY
# =============================================================================

def save_results(output: str,
                 processed: pd.DataFrame,
                 stats: Dict,
                 intervals: List[TimeInterval],
                 thresholds: Tuple[float, ...]):
    """Export all computed results to a structured Excel workbook."""
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # All data
        processed.to_excel(writer, sheet_name='C0_C7', index=False)

        # Interval statistics
        if 'interval_stats' in stats:
            stats['interval_stats'].to_excel(writer, sheet_name='tau_max_summary', index=False)
        if 'threshold_stats' in stats:
            stats['threshold_stats'].to_excel(writer, sheet_name='tau_max_thresholds', index=False)

        # Per-interval sheets and thresholded subsets
        for interval in intervals:
            mask = (
                (processed['Time'] >= pd.to_datetime(interval.start_date)) &
                (processed['Time'] <= pd.to_datetime(interval.end_date))
            )
            subset = processed.loc[mask]
            if not subset.empty:
                subset.to_excel(writer, sheet_name=interval.sheet_name[:31], index=False)
                for thr in thresholds:
                    sel = subset[subset['tau_max'] > thr]
                    if not sel.empty:
                        sheet_name = f"{interval.sheet_name}_tau_max>{thr:.3f}"[:31]
                        sel.to_excel(writer, sheet_name=sheet_name, index=False)


# =============================================================================
# 6. MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compute shear stresses, wave–current metrics, u_star, and δ_mix."
    )
    parser.add_argument('-i', '--input', default=CONFIG.input_file)
    parser.add_argument('-o', '--output', default=CONFIG.output_file)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    df = pd.read_excel(args.input)

    # Normalize possible tau_cr column names from Excel
    df = df.rename(columns={'Tau_cr': 'tau_cr', 'TAU_CR': 'tau_cr'})

    # Basic input integrity check
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required input columns: {missing}")

    processed = process_oceanographic_data(df)
    stats = calculate_interval_statistics(processed, TIME_INTERVALS, CONFIG.tau_max_thresholds)
    save_results(args.output, processed, stats, TIME_INTERVALS, CONFIG.tau_max_thresholds)

    print(f"Processing complete. Results saved to {args.output}")


if __name__ == "__main__":
    main()
