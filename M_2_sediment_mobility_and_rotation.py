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

import numpy as np
import pandas as pd
import logging
import sys
from pathlib import Path
from dataclasses import dataclass

import matplotlib.pyplot as plt  # for the rotation plot


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class MasterFileConfig:
    input_file: str = "INPUT_SEDIMENT_PARAMETERS.xlsx"
    output_file: str = "OUTPUT_SEDIMENT_MOBILITY_PARAMETERS.xlsx"
    rotation_plot: str = "COORDINATE_ROTATION_DIAGRAM.jpg"


@dataclass
class RotationConfig:
    x_col: str = "x"
    y_col: str = "y"
    rotated_x_col: str = "as"
    rotated_y_col: str = "cs"
    angle_degrees: float = -22.0   # clockwise rotation


@dataclass
class SedimentConfig:
    d50_col: str = "D50"
    rho_sed_col: str = "rho_sed"
    rho_sw_col: str = "rho_sw"
    nu_col: str = "nu"
    g_col: str = "g"
    sed_camp_col: str = "Sed_Camp"
    sample_col: str = "Sample_GIS"
    local_col: str = "Local"
    date_col: str = "date"


MASTER = MasterFileConfig()
ROT_CFG = RotationConfig()
SED_CFG = SedimentConfig()


# =============================================================================
# LOGGING
# =============================================================================

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )


# =============================================================================
# PHYSICAL FUNCTIONS
# =============================================================================

def rotate_coordinates(x, y, angle_deg, origin):
    """
    Rotate spatial coordinates (x, y) by angle_deg around origin = (x0, y0).

    angle_deg < 0 -> clockwise rotation.
    """
    theta = np.radians(angle_deg)
    x_t, y_t = x - origin[0], y - origin[1]
    x_r = x_t * np.cos(theta) - y_t * np.sin(theta)
    y_r = x_t * np.sin(theta) + y_t * np.cos(theta)
    return x_r + origin[0], y_r + origin[1]


def calculate_dimensionless_grain_size(d50, rho_s, rho_w, nu, g):
    """Soulsby (1997) dimensionless grain size D*."""
    return d50 * ((g * (rho_s - rho_w) / (rho_w * nu**2)) ** (1.0 / 3.0))


def calculate_critical_shields_parameter(D_star):
    """Soulsby (1997) critical Shields parameter θ_cr."""
    return 0.30 / (1.0 + 1.2 * D_star) + 0.055 * (1.0 - np.exp(-0.020 * D_star))


def calculate_critical_shear_stress(theta_cr, d50, rho_s, rho_w, g):
    """Critical shear stress τ_cr."""
    return theta_cr * (rho_s - rho_w) * g * d50


def calculate_critical_shear_velocity(theta_cr, d50, rho_s, rho_w, g):
    """Critical shear velocity u*_cr."""
    return np.sqrt(theta_cr * g * (rho_s - rho_w) * d50 / rho_w)


def calculate_settling_velocity(d50, rho_s, rho_w, nu, g):
    """Soulsby (1997) settling velocity w_s."""
    D_star = calculate_dimensionless_grain_size(d50, rho_s, rho_w, nu, g)
    return (nu / d50) * (np.sqrt(10.36**2 + 1.049 * D_star**3) - 10.36)


# =============================================================================
# PROCESSING
# =============================================================================

def process_sediment_sheet(df):
    """Compute D*, θ_cr, τ_cr, u*_cr, w_s for one campaign sheet."""
    df = df.copy()

    df["D_star"] = df.apply(
        lambda r: calculate_dimensionless_grain_size(
            r[SED_CFG.d50_col],
            r[SED_CFG.rho_sed_col],
            r[SED_CFG.rho_sw_col],
            r[SED_CFG.nu_col],
            r[SED_CFG.g_col],
        ),
        axis=1,
    )

    df["Theta_cr"] = df["D_star"].apply(calculate_critical_shields_parameter)

    df["Tau_cr"] = df.apply(
        lambda r: calculate_critical_shear_stress(
            r["Theta_cr"],
            r[SED_CFG.d50_col],
            r[SED_CFG.rho_sed_col],
            r[SED_CFG.rho_sw_col],
            r[SED_CFG.g_col],
        ),
        axis=1,
    )

    df["u_star_cr"] = df.apply(
        lambda r: calculate_critical_shear_velocity(
            r["Theta_cr"],
            r[SED_CFG.d50_col],
            r[SED_CFG.rho_sed_col],
            r[SED_CFG.rho_sw_col],
            r[SED_CFG.g_col],
        ),
        axis=1,
    )

    df["w_s"] = df.apply(
        lambda r: calculate_settling_velocity(
            r[SED_CFG.d50_col],
            r[SED_CFG.rho_sed_col],
            r[SED_CFG.rho_sw_col],
            r[SED_CFG.nu_col],
            r[SED_CFG.g_col],
        ),
        axis=1,
    )

    return df


# =============================================================================
# PLOT: COORDINATE ROTATION
# =============================================================================

def create_coordinate_rotation_plot(df, origin, output_path):
    """
    Create a QA plot showing original (x,y) and rotated (as,cs) coordinates,
    centered on the tracer origin (so Tracer ≈ (0,0) in the plot).
    """
    try:
        x0, y0 = origin

        if ROT_CFG.x_col not in df.columns or ROT_CFG.y_col not in df.columns:
            logging.warning("Cannot create rotation plot: x or y columns missing.")
            return

        if ROT_CFG.rotated_x_col not in df.columns or ROT_CFG.rotated_y_col not in df.columns:
            logging.warning("Cannot create rotation plot: rotated columns as/cs missing.")
            return

        x_orig = df[ROT_CFG.x_col] - x0
        y_orig = df[ROT_CFG.y_col] - y0

        x_rot = df[ROT_CFG.rotated_x_col] - x0
        y_rot = df[ROT_CFG.rotated_y_col] - y0

        fig, ax = plt.subplots(figsize=(7, 7))

        ax.scatter(x_orig, y_orig, alpha=0.6, s=35, label="Original (x, y)")
        ax.scatter(x_rot, y_rot, alpha=0.6, s=35, marker="x", label="Rotated (as, cs)")

        # Tracer origin at (0,0)
        ax.scatter(0.0, 0.0, marker="*", s=120, label="Tracer origin (injection)")

        ax.axhline(0, color="k", linewidth=0.6, alpha=0.5)
        ax.axvline(0, color="k", linewidth=0.6, alpha=0.5)

        ax.set_xlabel("Along-shore / cross-shore frame (m)\ncentered on Tracer")
        ax.set_ylabel("Along-shore / cross-shore frame (m)")
        ax.set_title("Coordinate rotation around Tracer (−22° clockwise)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal", adjustable="box")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Rotation plot saved to: {output_path}")

    except Exception as e:
        logging.error(f"Error creating rotation plot: {e}")


# =============================================================================
# STATISTICS
# =============================================================================

def calculate_summary_statistics(df):
    """Campaign and global stats for D50, Tau_cr, u_star_cr, w_s (and delta_mix if present).

    - Campaign stats are computed PER SITE (Local) and PER campaign (Sed_Camp), so Tavira/Vilamoura
      are not mixed together.
    - Also returns a 'campaign_samples' table with one row per sample.
    """
    stats = {}

    # A tidy per-sample table (what many users expect to see under 'campaign stats')
    keep_cols = [
        c for c in [
            SED_CFG.local_col,
            SED_CFG.sed_camp_col,
            SED_CFG.date_col,
            SED_CFG.sample_col,
            ROT_CFG.x_col,
            ROT_CFG.y_col,
            SED_CFG.d50_col,
            "Tau_cr",
            "u_star_cr",
            "w_s",
        ] if c in df.columns
    ]
    stats["campaign_samples"] = df[keep_cols].copy() if keep_cols else df.copy()

    # Grouped stats (mean/std/n)
    group_cols = []
    if SED_CFG.local_col in df.columns:
        group_cols.append(SED_CFG.local_col)
    if SED_CFG.sed_camp_col in df.columns:
        group_cols.append(SED_CFG.sed_camp_col)

    if group_cols:
        agg_map = {
            SED_CFG.d50_col: ["count", "mean", "std"],
            "Tau_cr": ["mean", "std"],
            "u_star_cr": ["mean", "std"],
            "w_s": ["mean", "std"],
        }
        stats["campaign_stats"] = df.groupby(group_cols).agg(agg_map)

        if "delta_mix" in df.columns:
            stats["max_delta_mix"] = (
                df.groupby(group_cols)["delta_mix"]
                .max()
                .reset_index()
            )

    stats["global_stats"] = {
        "total_samples": int(len(df)),
        "mean_D50": float(df[SED_CFG.d50_col].mean()) if SED_CFG.d50_col in df.columns else np.nan,
    }

    return stats


# =============================================================================
# SAVE RESULTS
# =============================================================================

def save_results(path, df, stats):
    """Write results to Excel workbook."""
    with pd.ExcelWriter(path) as writer:
        # 1) Full parameters table
        df.to_excel(
            excel_writer=writer,
            sheet_name="Sediment_Parameters",
            index=False,
        )

        # 2) Rotated coordinates sheet – only include columns that exist
        rot_cols = [
            col
            for col in [
                SED_CFG.sample_col,      # "Sample" (may be absent)
                SED_CFG.sed_camp_col,    # "Sed_Camp" (may be absent)
                ROT_CFG.x_col,           # "x"
                ROT_CFG.y_col,           # "y"
                ROT_CFG.rotated_x_col,   # "as"
                ROT_CFG.rotated_y_col,   # "cs"
            ]
            if col in df.columns
        ]

        if rot_cols:
            df[rot_cols].to_excel(
                excel_writer=writer,
                sheet_name="Rotated_Coordinates",
                index=False,
            )

        # 3) Campaign statistics (per-sample table)
        # Keep this sheet clean: one row per sample (no summary tables appended).
        if "campaign_samples" in stats:
            stats["campaign_samples"].to_excel(
                excel_writer=writer,
                sheet_name="Campaign_Statistics",
                index=False,
                startrow=0,
            )

        # 3b) Summary statistics (grouped mean/std) in a dedicated sheet
        if "campaign_stats" in stats:
            stats["campaign_stats"].to_excel(
                excel_writer=writer,
                sheet_name="Summary_Statistics",
            )

        # 4) Max mixing thickness (only if delta_mix present)
        if "max_delta_mix" in stats:
            stats["max_delta_mix"].to_excel(
                excel_writer=writer,
                sheet_name="Max_Mixing_Thickness",
                index=False,
            )

        # 5) Global statistics
        pd.DataFrame([stats["global_stats"]]).to_excel(
            excel_writer=writer,
            sheet_name="Global_Statistics",
            index=False,
        )


# =============================================================================
# MAIN
# =============================================================================

def main():
    setup_logging()
    logging.info("Processing sediment datasets...")

    # Read workbook
    xls = pd.ExcelFile(MASTER.input_file)

    # ------------------------------------------------------------------
    # Determine rotation origin(s) from Tracer sheet(s)
    # ------------------------------------------------------------------
    tracer_sheets = [s for s in xls.sheet_names if str(s).strip().lower().startswith("tracer")]

    if not tracer_sheets:
        logging.error('No sheet starting with "Tracer" was found in the input workbook. '
                      'Rotation origin cannot be defined.')
        sys.exit(1)

    origin_map = {}  # site_token -> (x0, y0)
    for sheet in tracer_sheets:
        try:
            tracer_df = pd.read_excel(xls, sheet_name=sheet)
        except Exception as e:
            logging.warning(f"Could not read tracer sheet '{sheet}': {e}")
            continue

        if ROT_CFG.x_col not in tracer_df.columns or ROT_CFG.y_col not in tracer_df.columns:
            logging.warning(
                f"Tracer sheet '{sheet}' is missing columns '{ROT_CFG.x_col}' and/or '{ROT_CFG.y_col}'. Skipping it."
            )
            continue

        # Site token rule:
        #   - If sheet name contains an underscore, use the part after the first underscore (e.g., Tracer_Vilamoura -> Vilamoura)
        #   - Otherwise use the full sheet name.
        site_token = sheet.split("_")[-1] if "_" in sheet else sheet
        origin_map[site_token] = (
            float(tracer_df[ROT_CFG.x_col].mean()),
            float(tracer_df[ROT_CFG.y_col].mean()),
        )
        logging.info(
            f"Rotation origin set from '{sheet}' at (x={origin_map[site_token][0]:.2f}, y={origin_map[site_token][1]:.2f})"
        )

    if not origin_map:
        logging.error('Tracer sheet(s) were found, but none contained valid "x" and "y" columns.')
        sys.exit(1)

    # Default origin if we cannot match a sheet to a site token
    default_origin = next(iter(origin_map.values()))
    default_site = next(iter(origin_map.keys()))

    processed_all = []
    processed_sheet_names = []
    skipped_sheet_info = []

    # Columns required in each sheet to perform sediment + rotation computations
    required_cols = [
        ROT_CFG.x_col,
        ROT_CFG.y_col,
        SED_CFG.d50_col,
        SED_CFG.rho_sed_col,
        SED_CFG.rho_sw_col,
        SED_CFG.nu_col,
        SED_CFG.g_col,
    ]

    for sheet in xls.sheet_names:
        try:
            df = pd.read_excel(xls, sheet_name=sheet)

            # Check required columns
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                msg = f"Sheet '{sheet}' skipped: missing required columns {missing}"
                logging.warning(msg)
                skipped_sheet_info.append(msg)
                continue

            # If Sed_Camp does not exist, create it from sheet name
            if SED_CFG.sed_camp_col not in df.columns:
                df[SED_CFG.sed_camp_col] = sheet

            # Add a 'Site' column (helps later grouping/plotting)
            # Prefer an existing 'Local' column; otherwise infer from sheet name token.
            if "Site" not in df.columns:
                if "Local" in df.columns:
                    df["Site"] = df["Local"].astype(str)
                else:
                    df["Site"] = sheet.split("_")[-1] if "_" in sheet else sheet

            # Drop rows where *all* required numeric fields are NaN
            df_non_empty = df.dropna(subset=required_cols, how="all")
            if df_non_empty.empty:
                msg = f"Sheet '{sheet}' skipped: all rows have NaN in required fields."
                logging.warning(msg)
                skipped_sheet_info.append(msg)
                continue

            # Enforce expected campaign based on sheet name to avoid mixing duplicated rows
            # (e.g., a Tracer_* sheet accidentally containing Av_Native_Sed rows).
            sheet_lower = sheet.lower()
            if sheet_lower.startswith("tracer"):
                expected_camp = "Tracer"
            elif sheet_lower.startswith("nat_sed") or sheet_lower.startswith("native"):
                expected_camp = "Av_Native_Sed"
            else:
                expected_camp = None

            if expected_camp and SED_CFG.sed_camp_col in df_non_empty.columns:
                df_non_empty = df_non_empty[df_non_empty[SED_CFG.sed_camp_col] == expected_camp].copy()

            if expected_camp and SED_CFG.sed_camp_col not in df_non_empty.columns:
                df_non_empty[SED_CFG.sed_camp_col] = expected_camp

            # If the sheet got emptied by filtering, skip it
            if df_non_empty.empty:
                msg = f"Sheet '{sheet}' skipped: no rows match expected campaign '{expected_camp}'."
                logging.warning(msg)
                skipped_sheet_info.append(msg)
                continue

            df_proc = process_sediment_sheet(df_non_empty)

            # Decide which origin to use for this sheet:
            # If the sheet name is like Something_Vilamoura, match 'Vilamoura' to a tracer token.
            site_token = sheet.split("_")[-1] if "_" in sheet else sheet
            origin = origin_map.get(site_token, default_origin)
            origin_site = site_token if site_token in origin_map else default_site

            # Ensure Local matches the site token from the sheet name (helps avoid mixed group stats)
            df_proc[SED_CFG.local_col] = origin_site

            df_proc["Origin_Site"] = origin_site
            df_proc["Origin_x"] = origin[0]
            df_proc["Origin_y"] = origin[1]

            # Rotate coordinates (x,y → as, cs) w.r.t. site-specific Tracer origin
            df_proc[[ROT_CFG.rotated_x_col, ROT_CFG.rotated_y_col]] = df_proc.apply(
                lambda r: rotate_coordinates(
                    r[ROT_CFG.x_col],
                    r[ROT_CFG.y_col],
                    ROT_CFG.angle_degrees,
                    origin,
                ),
                axis=1,
            ).tolist()

            processed_all.append(df_proc)
            processed_sheet_names.append(sheet)
            logging.info(f"Sheet '{sheet}': processed {len(df_proc)} valid samples (origin: {origin_site}).")

        except Exception as e:
            msg = f"Sheet '{sheet}' skipped due to error: {e}"
            logging.error(msg)
            skipped_sheet_info.append(msg)
            continue

    if not processed_all:
        logging.error("No valid sheets were processed. Check required columns and data content.")
        if skipped_sheet_info:
            logging.error("Skipped sheets summary:")
            for msg in skipped_sheet_info:
                logging.error("  - " + msg)
        sys.exit(1)

    final_df = pd.concat(processed_all, ignore_index=True)

    logging.info("SUMMARY OF SHEET PROCESSING:")
    logging.info(f"  Sheets processed ({len(processed_sheet_names)}): {processed_sheet_names}")
    if skipped_sheet_info:
        logging.info("  Sheets skipped:")
        for msg in skipped_sheet_info:
            logging.info("    - " + msg)
    logging.info(f"  Total valid samples in final dataset: {len(final_df)}")

    # ------------------------------------------------------------------
    # Create rotation QA plot(s): one per origin/site when possible
    # ------------------------------------------------------------------
    try:
        # If the user has multiple tracer origins, produce one plot per origin site token.
        for site, origin in origin_map.items():
            df_site = final_df[final_df["Origin_Site"] == site] if "Origin_Site" in final_df.columns else final_df
            if df_site.empty:
                continue

            plot_name = Path(MASTER.rotation_plot)
            plot_path = (
                plot_name.with_name(f"{plot_name.stem}_{site}{plot_name.suffix}")
                if len(origin_map) > 1
                else plot_name
            )

            create_coordinate_rotation_plot(
                df=df_site,
                origin=origin,
                output_path=str(plot_path),
            )
    except Exception as e:
        logging.warning(f"Rotation plot step skipped due to error: {e}")

    # ------------------------------------------------------------------
    # Stats + Excel output
    # ------------------------------------------------------------------
    stats = calculate_summary_statistics(final_df)
    save_results(Path(MASTER.output_file), final_df, stats)

    logging.info("✔ Completed successfully.")


if __name__ == "__main__":
    main()
