"""
===================================
Directional statistics script
File:  tracer_metrics_processing.py
===================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: January 2026

Computes exceedance-weighted and unweighted directional metrics used to
compare hydrodynamic forcing with observed tracer centre-of-mass (TCOM)
displacement.

Paper-only behaviour
--------------------
- Retains all calculations required to derive exceedance-driven directions
  robustly (integrated exceedance intensity, mean resultant lengths,
  sine/cosine components, exceedance-weighted circular means, and merged
  exceedance directions).
- Exports ONLY the manuscript-facing alignment table (one row per sampling
  interval), containing:

    Sampling intervals
    TCOMd (m)
    TCOMdir (deg)
    CdirAll (deg)
    CdirEm (deg)
    Diff. (TCOMdir:CdirAll) (deg)
    Diff. (TCOMdir:CdirEm) (deg)

Threshold selection for the reported exceedance direction (CdirEm)
-----------------------------------------------------------------
For each sampling interval, the script selects a single critical shear-stress
threshold (tau_cr) as follows:
  1) Prefer PREFERRED_TAU_CR if CdirEm is defined (i.e. exceedance occurs and
     yields a meaningful direction).
  2) Otherwise, fall back to FALLBACK_TAU_CR.

Directional conventions
-----------------------
- Oceanographic azimuth convention: 0° = North, 90° = East.
- Minimum angular differences are reported in the range [0°, 180°].
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Union

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG (edit paths/intervals here)
# =============================================================================
HERE = Path(__file__).parent.resolve()

INPUT_WORKBOOK = HERE / "INPUT_using_OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx"
SHEET_NAME = "C0_C7"  # must contain columns: Time, ve, vn, tau_max

TRACER_CM_PATH = HERE / "TCOM_coordinates.xlsx"  # must contain COM labels + x,y

OUTPUT_WORKBOOK = HERE / "OUTPUT_TCOM_ALIGNMENT_MANUSCRIPT.xlsx"
OUTPUT_SHEET = "TCOM_alignment_table_final"

# Shear-stress thresholds and exceedance exponents
TAU_CR_VALUES = (0.259, 0.303)  # N/m^2
N_VALUES = (1.0, 2.0)

# What the paper reports (preferred then fallback)
PREFERRED_TAU_CR = 0.303
FALLBACK_TAU_CR = 0.259

# Sampling intervals
@dataclass
class TimeInterval:
    start_date: str
    end_date: str
    name: str  # e.g. 'C0_C2'


TIME_INTERVALS: List[TimeInterval] = [
    TimeInterval("2020-08-18 17:50:08", "2020-08-30 23:50:08", "C0_C2"),
    TimeInterval("2020-08-31 00:05:08", "2020-09-28 23:50:08", "C2_C3"),
    TimeInterval("2020-09-29 00:05:08", "2020-11-14 23:50:08", "C3_C4"),
    TimeInterval("2020-11-15 00:05:08", "2021-01-14 23:50:08", "C4_C5"),
    TimeInterval("2021-01-15 00:05:08", "2021-04-22 23:50:08", "C5_C6"),
    TimeInterval("2021-04-23 00:05:08", "2021-10-27 23:50:08", "C6_C7"),
]


# =============================================================================
# Helpers: time handling and angles
# =============================================================================

def _ensure_time(df: pd.DataFrame) -> pd.DataFrame:
    if "Time" not in df.columns:
        raise ValueError("Input must include a 'Time' column.")
    if not np.issubdtype(df["Time"].dtype, np.datetime64):
        df = df.copy()
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
    return df


def _forward_dt_seconds(df: pd.DataFrame) -> pd.Series:
    df = _ensure_time(df.sort_values("Time"))
    dt_fwd = (df["Time"].shift(-1) - df["Time"]).dt.total_seconds()
    dt_bwd = (df["Time"] - df["Time"].shift(1)).dt.total_seconds()
    dt = dt_fwd.fillna(dt_bwd).clip(lower=0.0)
    return dt


def _azimuth_from_components_scalar(ve: float, vn: float) -> float:
    """Oceanographic azimuth for scalar components (0=N, 90=E)."""
    ang = math.degrees(math.atan2(ve, vn))
    return (ang + 360.0) % 360.0


def _azimuth_from_components_vec(ve: np.ndarray, vn: np.ndarray) -> np.ndarray:
    """Oceanographic azimuth for vector components (0=N, 90=E)."""
    ang = np.degrees(np.arctan2(ve, vn))
    return (ang + 360.0) % 360.0


def _min_angle_diff_deg(a_deg: float, b_deg: float) -> float:
    """Smallest absolute angular separation between two azimuths (0..360)."""
    if a_deg is None or b_deg is None or pd.isna(a_deg) or pd.isna(b_deg):
        return np.nan
    d = float((a_deg - b_deg) % 360.0)
    return 360.0 - d if d > 180.0 else d


def subset_interval(df: pd.DataFrame, start: str, end: str) -> pd.DataFrame:
    df = _ensure_time(df)
    start_dt = pd.to_datetime(start)
    end_dt = pd.to_datetime(end)
    return df[(df["Time"] >= start_dt) & (df["Time"] <= end_dt)].copy()


# =============================================================================
# Core directional computations
# =============================================================================

def all_times_component_mean(df: pd.DataFrame) -> dict:
    """Background direction (CdirAll) from dt-weighted component means."""
    if df.empty:
        return {"dir_deg": np.nan, "mag_mps": np.nan}

    df = _ensure_time(df)
    dt = _forward_dt_seconds(df)

    ve = df["ve"].to_numpy(float)
    vn = df["vn"].to_numpy(float)

    w = np.where(np.isfinite(dt), dt.to_numpy(float), 0.0)
    wsum = float(np.nansum(w))
    if wsum <= 0.0:
        return {"dir_deg": np.nan, "mag_mps": np.nan}

    Ue = float(np.nansum(ve * w) / wsum)
    Un = float(np.nansum(vn * w) / wsum)
    mag = float(np.hypot(Ue, Un))

    if np.isclose(mag, 0.0, atol=1e-15):
        return {"dir_deg": np.nan, "mag_mps": 0.0}

    return {"dir_deg": _azimuth_from_components_scalar(Ue, Un), "mag_mps": mag}


def exceedance_weighted_direction(df: pd.DataFrame, tau_cr: float, n: float) -> dict:
    """Exceedance-weighted circular mean direction and coherence.

    Returns:
      theta_bar_deg : exceedance-weighted mean direction
      R             : mean resultant length (0..1)
      S_E_sum/S_N_sum: weighted sine/cosine sums
      IE_n          : integrated exceedance intensity (sum of weights)
      N             : number of exceedance samples
    """
    if df.empty:
        return {
            "theta_bar_deg": np.nan,
            "R": 0.0,
            "S_E_sum": 0.0,
            "S_N_sum": 0.0,
            "IE_n": 0.0,
            "N": 0,
        }

    df = _ensure_time(df).sort_values("Time")
    dt = _forward_dt_seconds(df).to_numpy(float)

    exceed_mask = (df["tau_max"].to_numpy(float) > float(tau_cr))
    if not exceed_mask.any():
        return {
            "theta_bar_deg": np.nan,
            "R": 0.0,
            "S_E_sum": 0.0,
            "S_N_sum": 0.0,
            "IE_n": 0.0,
            "N": 0,
        }

    tau = df["tau_max"].to_numpy(float)
    exceed = np.clip(tau - float(tau_cr), a_min=0.0, a_max=None)

    # weights: [(tau_max - tau_cr)_+]^n * dt
    w = (exceed ** float(n)) * dt

    # keep only exceedance samples
    w = w[exceed_mask]
    ve = df.loc[exceed_mask, "ve"].to_numpy(float)
    vn = df.loc[exceed_mask, "vn"].to_numpy(float)

    if not np.isfinite(w).any() or float(np.nansum(w)) <= 0.0:
        return {
            "theta_bar_deg": np.nan,
            "R": 0.0,
            "S_E_sum": 0.0,
            "S_N_sum": 0.0,
            "IE_n": 0.0,
            "N": int(exceed_mask.sum()),
        }

    theta_deg = _azimuth_from_components_vec(ve, vn)
    theta_rad = np.radians(theta_deg)

    # unit vectors on circle (east component = sin, north component = cos)
    e_unit = np.sin(theta_rad)
    n_unit = np.cos(theta_rad)

    S_E = float(np.nansum(w * e_unit))
    S_N = float(np.nansum(w * n_unit))
    IE_n = float(np.nansum(w))

    if IE_n <= 0.0 or (np.isclose(S_E, 0.0, atol=1e-15) and np.isclose(S_N, 0.0, atol=1e-15)):
        return {
            "theta_bar_deg": np.nan,
            "R": 0.0,
            "S_E_sum": S_E,
            "S_N_sum": S_N,
            "IE_n": IE_n,
            "N": int(exceed_mask.sum()),
        }

    theta_bar = float((math.degrees(math.atan2(S_E, S_N)) + 360.0) % 360.0)
    R = float(np.hypot(S_E, S_N) / IE_n)

    return {
        "theta_bar_deg": theta_bar,
        "R": R,
        "S_E_sum": S_E,
        "S_N_sum": S_N,
        "IE_n": IE_n,
        "N": int(exceed_mask.sum()),
    }


def _circ_mean_two_weighted(dir1_deg, w1, dir2_deg, w2) -> float:
    """Weighted circular mean of two azimuths (degrees)."""
    angles = []
    weights = []

    if dir1_deg is not None and not pd.isna(dir1_deg) and w1 is not None and np.isfinite(w1) and w1 > 0:
        angles.append(np.deg2rad(float(dir1_deg)))
        weights.append(float(w1))
    if dir2_deg is not None and not pd.isna(dir2_deg) and w2 is not None and np.isfinite(w2) and w2 > 0:
        angles.append(np.deg2rad(float(dir2_deg)))
        weights.append(float(w2))

    if not angles:
        return np.nan
    if len(angles) == 1:
        return float((np.rad2deg(angles[0]) + 360.0) % 360.0)

    a = np.asarray(angles, float)
    w = np.asarray(weights, float)

    E = float(np.sum(w * np.sin(a)))
    N = float(np.sum(w * np.cos(a)))
    if np.isclose(E, 0.0, atol=1e-15) and np.isclose(N, 0.0, atol=1e-15):
        return np.nan

    return _azimuth_from_components_scalar(E, N)


def merge_exceedance_directions(theta_n1: float, R_n1: float, IE_n1: float,
                                theta_n2: float, R_n2: float, IE_n2: float) -> dict:
    """Merged exceedance directions.

    - CdirEm uses R-weighting (directional coherence weighting).
    - Jdir (not exported) uses IE-weighting (energetic weighting).
    """
    theta_CdirEm = _circ_mean_two_weighted(theta_n1, R_n1, theta_n2, R_n2)
    theta_Jdir = _circ_mean_two_weighted(theta_n1, IE_n1, theta_n2, IE_n2)
    return {"CdirEm_deg": theta_CdirEm, "Jdir_deg": theta_Jdir}


# =============================================================================
# Tracer COM displacement (TCOMd, TCOMdir)
# =============================================================================

def make_tracer_disp_from_cm_excel(cm_path: Union[str, Path],
                                   intervals: Sequence[TimeInterval]) -> pd.DataFrame:
    cm = pd.read_excel(cm_path)
    lower = {c.lower(): c for c in cm.columns}

    cmass_col = lower.get("cmass")
    if cmass_col is None:
        cmass_col = next((c for c in cm.columns if "cm" in c.lower()), None)
    if cmass_col is None:
        raise ValueError("Could not find a COM label column (e.g. 'Cmass').")

    x_col = lower.get("x")
    if x_col is None:
        x_col = next((c for c in cm.columns if re.fullmatch(r"x\w*", c.lower())), None)
    y_col = lower.get("y")
    if y_col is None:
        y_col = next((c for c in cm.columns if re.fullmatch(r"y\w*", c.lower())), None)
    if x_col is None or y_col is None:
        raise ValueError("Could not find 'x'/'y' COM coordinate columns.")

    # map labels like C0, C2, ... to (x,y)
    mapping = {}
    for _, r in cm.iterrows():
        label = str(r[cmass_col])
        m = re.search(r"(\d+)", label)
        if not m:
            continue
        key = f"C{m.group(1)}"
        mapping[key] = (float(r[x_col]), float(r[y_col]))

    rows = []
    for ti in intervals:
        if "_" not in ti.name:
            rows.append({"Interval": ti.name, "TCOMd_m": np.nan, "TCOMdir_deg": np.nan})
            continue

        a, b = ti.name.split("_", 1)
        if a in mapping and b in mapping:
            x0, y0 = mapping[a]
            x1, y1 = mapping[b]
            dx = x1 - x0
            dy = y1 - y0
            dist = float(math.hypot(dx, dy))
            az = float((math.degrees(math.atan2(dx, dy)) + 360.0) % 360.0)  # 0=N, 90=E
            rows.append({"Interval": ti.name, "TCOMd_m": dist, "TCOMdir_deg": az})
        else:
            rows.append({"Interval": ti.name, "TCOMd_m": np.nan, "TCOMdir_deg": np.nan})

    return pd.DataFrame(rows)


# =============================================================================
# Per-interval diagnostics (computed for BOTH tau_cr values)
# =============================================================================

def compute_interval_diagnostics(processed_data: pd.DataFrame,
                                 time_intervals: Sequence[TimeInterval],
                                 tau_cr_values: Sequence[float],
                                 n_values: Sequence[float]) -> pd.DataFrame:
    records = []

    for ti in time_intervals:
        df_i = subset_interval(processed_data, ti.start_date, ti.end_date)
        logging.info("Interval %s: %d records", ti.name, len(df_i))

        cdir_all = all_times_component_mean(df_i)["dir_deg"]

        for tau_cr in tau_cr_values:
            ew1 = exceedance_weighted_direction(df_i, tau_cr=float(tau_cr), n=float(n_values[0]))
            ew2 = exceedance_weighted_direction(df_i, tau_cr=float(tau_cr), n=float(n_values[1]))

            merged = merge_exceedance_directions(
                theta_n1=ew1["theta_bar_deg"], R_n1=ew1["R"], IE_n1=ew1["IE_n"],
                theta_n2=ew2["theta_bar_deg"], R_n2=ew2["R"], IE_n2=ew2["IE_n"],
            )

            records.append({
                "Interval": ti.name,
                "tau_cr": float(tau_cr),
                "CdirAll_deg": cdir_all,
                "CdirEm_deg": merged["CdirEm_deg"],
                # internal diagnostics (not exported, but computed)
                "IE_n1": ew1["IE_n"],
                "IE_n2": ew2["IE_n"],
                "R_n1": ew1["R"],
                "R_n2": ew2["R"],
            })

    return pd.DataFrame.from_records(records)


def build_manuscript_alignment_table(merged_df: pd.DataFrame,
                                     tracer_df: pd.DataFrame,
                                     preferred_tau_cr: float,
                                     fallback_tau_cr: float,
                                     decimals: int = 2) -> pd.DataFrame:
    """Build the paper-facing table (one row per interval).

    IMPORTANT: Avoid using `or` with pandas Series (ambiguous truth). We
    explicitly check for None.
    """
    df = merged_df.merge(tracer_df, on="Interval", how="left")

    out_rows = []
    for interval, g in df.groupby("Interval", sort=False):
        pref = g.loc[g["tau_cr"].eq(float(preferred_tau_cr))]
        fb = g.loc[g["tau_cr"].eq(float(fallback_tau_cr))]

        def pick(block: pd.DataFrame) -> Union[pd.Series, None]:
            if block is None or block.empty:
                return None
            good = block.loc[block["CdirEm_deg"].notna()]
            if good.empty:
                return None
            return good.iloc[0]

        row = pick(pref)
        if row is None:
            row = pick(fb)
        if row is None:
            row = g.iloc[0]

        tcomd = row.get("TCOMd_m", np.nan)
        tcomdir = row.get("TCOMdir_deg", np.nan)
        cdirall = row.get("CdirAll_deg", np.nan)
        cdirem = row.get("CdirEm_deg", np.nan)

        out_rows.append({
            "Sampling intervals": str(interval).replace("_", "–"),
            "TCOMd (m)": tcomd,
            "TCOMdir (º)": tcomdir,
            "CdirAll (º)": cdirall,
            "CdirEm (º)": cdirem,
            "Diff. (TCOMdir:CdirAll) (º)": _min_angle_diff_deg(tcomdir, cdirall),
            "Diff. (TCOMdir:CdirEm) (º)": _min_angle_diff_deg(tcomdir, cdirem),
        })

    out = pd.DataFrame(out_rows)

    # Round numeric columns for manuscript presentation
    for c in out.columns:
        if c == "Sampling intervals":
            continue
        out[c] = pd.to_numeric(out[c], errors="coerce").round(decimals)

    return out


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not INPUT_WORKBOOK.exists():
        raise FileNotFoundError(f"Missing input workbook: {INPUT_WORKBOOK}")

    xls = pd.ExcelFile(INPUT_WORKBOOK)
    if SHEET_NAME not in xls.sheet_names:
        raise ValueError(
            f"Sheet '{SHEET_NAME}' not found in {INPUT_WORKBOOK.name}. "
            f"Available: {xls.sheet_names}"
        )

    data = pd.read_excel(xls, sheet_name=SHEET_NAME, parse_dates=["Time"])
    required = {"Time", "ve", "vn", "tau_max"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Input data missing required columns: {missing}")

    if not TRACER_CM_PATH.exists():
        raise FileNotFoundError(f"Missing tracer COM workbook: {TRACER_CM_PATH}")

    tracer_df = make_tracer_disp_from_cm_excel(TRACER_CM_PATH, TIME_INTERVALS)
    merged_df = compute_interval_diagnostics(data, TIME_INTERVALS, TAU_CR_VALUES, N_VALUES)

    table = build_manuscript_alignment_table(
        merged_df=merged_df,
        tracer_df=tracer_df,
        preferred_tau_cr=PREFERRED_TAU_CR,
        fallback_tau_cr=FALLBACK_TAU_CR,
        decimals=2,
    )

    with pd.ExcelWriter(OUTPUT_WORKBOOK, engine="openpyxl", mode="w") as writer:
        table.to_excel(writer, sheet_name=OUTPUT_SHEET, index=False)

    logging.info("Wrote manuscript table: %s (sheet: %s)", OUTPUT_WORKBOOK.name, OUTPUT_SHEET)


if __name__ == "__main__":
    main()
