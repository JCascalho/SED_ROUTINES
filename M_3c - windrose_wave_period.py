"""
===============================================
Windrose plot generator for wave period data
File: windrose_wave_period.py
===============================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Description
-----------
This script generates a windrose-style plot showing the directional
distribution of wave period, using either the zero-crossing period (Tz)
or the peak period (Tp). Directions are binned into 32 canonical compass
sectors (11.25° resolution). Output figures are publication-grade and can
be exported in vector (PDF, SVG) or raster (PNG, JPG) formats at high
resolution.

Frequencies can be normalized either globally (relative to all
observations) or per direction, such that the frequency classes within
each directional sector sum to 100%. The resulting visualization provides
a compact overview of the directional structure of wave periods relevant
to wave–current interactions and sediment dynamics.

Input
-----
- Excel file (default):
    INPUT_using_OUTPUT_CURRENTS_WAVES_DATA_BASE_WITH_DELTA_MIX.xlsx
- Sheet name (command-line argument --sheet / -s, default: "C0_C7")
- Wave-period parameter (command-line argument --parameter / -p):
    - "Tz" : zero-crossing wave period [s]
    - "Tp" : peak wave period [s]
- Required columns:
    - dir_w : wave direction in degrees (0–360; 0° = North; “from” convention)
    - Tz or Tp (depending on the selected parameter)

Output
------
- Windrose plot file:
    (argument --output / -o, default: WINDROSE_Tz_C0_C7.jpg)
- Legend file:
    (argument --legend / -l, default: legend_Tz.jpg)
- Output formats controlled by:
    --format, --legend-format  (pdf, png, svg, jpg)

Usage examples
--------------
Basic usage (output format inferred from filename; Tz analysis):
  python S_7_6 - wind rose plot generator for wave period data.py \
      --sheet "C0_C7" --parameter Tz

Explicit JPG output with absolute paths:
  python S_7_6 - wind rose plot generator for wave period data.py \
      -s "C0_C7" --parameter Tz \
      -o "C:/Temp/WINDROSE_Tz.jpg" --format jpg \
      --legend "C:/Temp/legend_Tz.jpg" --legend-format jpg

Per-direction normalization (each sector sums to 100%; Tp analysis):
  python S_7_6 - wind rose plot generator for wave period data.py \
      -s "C0_C7" --parameter Tp \
      -o "WINDROSE_Tp.jpg" --format jpg --normalize per-direction
"""


import argparse
import os
import sys
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager

# ==================== CONFIGURATION ====================

INPUT_FILE = "INPUT_using_OUTPUT_CURRENTS_WAVES_DATA_BASE.xlsx"

# Wave period bins (seconds) and class names/colors
WAVE_PERIOD_BINS = [0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, np.inf]
WAVE_PERIOD_CLASS = [
    "<3 s",
    "3–6 s",
    "6–9 s",
    "9–12 s",
    "12–15 s",
    "15–18 s",
    "≥18 s",
]
WAVE_PERIOD_COLORS = [
    (0/255, 166/255, 255/255),  # Very short - light blue
    (0/255, 183/255, 166/255),  # Short - teal
    (255/255, 242/255, 0/255),  # Moderate - yellow
    (255/255, 161/255, 0/255),  # Long - orange
    (255/255, 0/255, 0/255),    # Very long - red
    (153/255, 102/255, 51/255), # Swell dominated - brown
    (51/255, 51/255, 51/255),   # Extreme swell - dark gray
]

# 32 directional sectors (11.25° each) – canonical 32-point compass labels
DIRECTION_BINS = np.linspace(0, 360, 33)
DIRECTION_LABELS = [
    "N", "NbE", "NNE", "NEbN", "NE", "NEbE", "ENE", "EbN",
    "E", "EbS", "ESE", "SEbE", "SE", "SEbS", "SSE", "SbE",
    "S", "SbW", "SSW", "SWbS", "SW", "SWbW", "WSW", "WbS",
    "W", "WbN", "WNW", "NWbW", "NW", "NWbN", "NNW", "NbW",
]

PLOT_CONFIG = {
    "figsize": (12, 12),
    "dpi": 600,  # for raster formats
    "grid_style": {"linestyle": "-", "color": "k", "linewidth": 0.8, "alpha": 0.35},
    "bar_style": {"edgecolor": "black", "linewidth": 0.8, "alpha": 1.0},
    "radial_percentages": [10, 20, 30, 40],  # reference rings (%)
    "direction_label_fontsize": 16,
    "direction_label_fontweight": "bold",
    "percentage_label_fontsize": 14,
    "title_fontsize": 18,
    "percent_label_angle_deg": 85,  # angle where % labels are placed
}

LEGEND_CONFIG = {
    "figsize": (6, 3.5),
    "font_properties": {"weight": "bold", "size": 12},
    "title": "Wave Period Classes",
    "edge_lw": 0.8,
}

# ==================== CLI ====================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality windrose plot for wave period (Tz/Tp).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --sheet "C0_C7"
  %(prog)s -s "C0_C7" -o results/WIND_ROSE_Tz.jpg --format jpg --normalize per-direction
""",
    )
    parser.add_argument("--sheet", "-s", default="C0_C7", #select the desired Excel file sheet 
                        help="Excel sheet name containing wave data (default: C0_C7)") #select the desired Excel file sheet 
    # Default to JPG filenames so Spyder/no-args users get JPGs by default
    parser.add_argument("--output", "-o", default="WIND_ROSE_Tz_C0_C7.jpg",
                        help="Output filename for windrose plot (default: WIND_ROSE_Tz_C0_C7.jpg)")
    parser.add_argument("--legend", "-l", default="legend_Tz.jpg",
                        help="Output filename for legend (default: legend_Tz.jpg)")
    parser.add_argument("--input", "-i", default=INPUT_FILE,
                        help=f"Input Excel file path (default: {INPUT_FILE})")
    parser.add_argument("--parameter", "-p", choices=["Tz", "Tp"], default="Tz",
                        help="Wave period parameter to analyze: 'Tz' (zero-crossing) or 'Tp' (peak).")
    parser.add_argument("--format", choices=["pdf", "png", "svg", "jpg"], default=None,
                        help="Override plot format (inferred from filename if omitted)")
    parser.add_argument("--legend-format", choices=["pdf", "png", "svg", "jpg"], default=None,
                        help="Override legend format (inferred from filename if omitted)")
    parser.add_argument("--normalize", choices=["global", "per-direction"], default="global",
                        help="Normalization: 'global' (default) or 'per-direction'")
    return parser.parse_args()

# ==================== SAVEFIG HELPERS (robust JPG support) ====================

def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)

def _post_save_report(path: str) -> None:
    apath = os.path.abspath(path)
    print(f"[savefig] attempted: {apath}")
    if os.path.exists(apath):
        print(f"[savefig] OK: file exists at {apath}")
    else:
        print(f"[savefig] WARNING: file not found at {apath}")

def _print_supported_filetypes():
    ft = plt.gcf().canvas.get_supported_filetypes()
    kinds = ", ".join(sorted(ft.keys()))
    print(f"[matplotlib] supported output formats: {kinds}")

def _safe_savefig(path: str, fig_format: str, dpi: int, bbox_inches: str = "tight") -> None:
    """
    Save a figure robustly:
    - Creates parent directories if needed.
    - Normalizes jpg/jpeg spelling.
    - If Pillow is missing, falls back to PNG automatically.
    """
    fmt = (fig_format or "").lower()
    fmt = "jpeg" if fmt in {"jpg", "jpeg"} else fmt
    _ensure_parent_dir(path)

    try:
        plt.savefig(path, dpi=dpi, bbox_inches=bbox_inches, format=fmt if fmt else None)
    except ValueError as e:
        # If jpg is requested but Pillow is missing, fallback to PNG
        if fmt == "jpeg":
            fallback_path = os.path.splitext(path)[0] + ".png"
            print(f"Warning: JPEG export requires Pillow. Falling back to PNG: {fallback_path}", file=sys.stderr)
            plt.savefig(fallback_path, dpi=dpi, bbox_inches=bbox_inches, format="png")
            _post_save_report(fallback_path)
            return
        raise
    _post_save_report(path)

def _figure_format(filename: str, arg_format: str, default: str = "jpg") -> str:
    # Default to 'jpg' so even with no args we try to write JPG
    ext = os.path.splitext(filename)[1].lower().lstrip(".")
    if arg_format:
        return arg_format
    if ext in {"pdf", "png", "svg", "jpg", "jpeg"}:
        return "jpg" if ext == "jpeg" else ext
    return default

# ==================== I/O & VALIDATION ====================

def validate_input_file(file_path: str) -> None:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"No read permission for file: {file_path}")

def load_wave_data(file_path: str, sheet_name: str, parameter_name: str) -> pd.DataFrame:
    xl = pd.ExcelFile(file_path)
    if sheet_name not in xl.sheet_names:
        raise ValueError(f"Sheet '{sheet_name}' not found. Available sheets: {xl.sheet_names}")
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    required = ["dir_w", parameter_name]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Available: {list(df.columns)}")

    # Basic checks
    if (df["dir_w"] < 0).any() or (df["dir_w"] > 360).any():
        print("Warning: Some 'dir_w' values outside 0–360°.", file=sys.stderr)
    if (df[parameter_name] < 0).any():
        raise ValueError(f"Negative {parameter_name} values found.")

    # Optional realism check (esp. for Tz)
    unrealistic = df[parameter_name] > 30
    if unrealistic.any():
        print(f"Warning: {unrealistic.sum()} {parameter_name} values > 30 s detected.", file=sys.stderr)

    # Drop NaNs
    start = len(df)
    df = df.dropna(subset=required)
    if df.empty:
        raise ValueError("No valid data remain after cleaning.")
    dropped = start - len(df)
    if dropped > 0:
        print(f"Note: dropped {dropped} rows with missing {required}.", file=sys.stderr)

    return df

# ==================== CORE LOGIC ====================

def process_wave_period_data(
    wave_directions: np.ndarray,
    wave_periods: np.ndarray,
    normalize: str = "global"
) -> Dict[str, List[float]]:
    """Bin directions and wave period classes, then normalize counts."""
    # Wrap to [0,360) and bin robustly; ensure 0° -> first sector
    dir_idx = np.digitize(wave_directions % 360.0, DIRECTION_BINS, right=False) - 1
    dir_idx = dir_idx % 32  # handle -1 and 32 safely

    per_idx = np.digitize(wave_periods, WAVE_PERIOD_BINS, right=False) - 1
    valid = (per_idx >= 0) & (per_idx < len(WAVE_PERIOD_BINS) - 1)

    counts = {lab: [0] * (len(WAVE_PERIOD_BINS) - 1) for lab in DIRECTION_LABELS}
    for d, p in zip(dir_idx[valid], per_idx[valid]):
        counts[DIRECTION_LABELS[d]][p] += 1

    # Normalize
    if normalize == "global":
        total = float(sum(sum(v) for v in counts.values()))
        if total <= 0:
            raise ValueError("No valid measurements found for processing.")
        for k in counts:
            counts[k] = [c / total for c in counts[k]]
    else:
        for k, vals in counts.items():
            s = float(sum(vals))
            counts[k] = [(v / s) if s > 0 else 0.0 for v in vals]

    return counts

# ==================== PLOTTING ====================

def create_windrose_plot(
    freq_data: Dict[str, List[float]],
    output_filename: str,
    fig_format: str,
    parameter_name: str,
    normalize: str
) -> None:
    fig = plt.figure(figsize=PLOT_CONFIG["figsize"])
    ax = fig.add_subplot(111, polar=True)

    # Oceanographic orientation: 0° at North, clockwise
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    theta_positions = np.linspace(0.0, 2 * np.pi, len(DIRECTION_LABELS), endpoint=False)
    sector_width = 2 * np.pi / len(DIRECTION_LABELS)

    # Grid
    ax.grid(True, **PLOT_CONFIG["grid_style"])

    # Stacked bars by period class
    for i, direction in enumerate(DIRECTION_LABELS):
        bottom = 0.0
        for cls in range(len(WAVE_PERIOD_CLASS)):
            frac = freq_data[direction][cls]
            if frac <= 0:
                continue
            ax.bar(theta_positions[i], frac, width=sector_width, bottom=bottom,
                   color=WAVE_PERIOD_COLORS[cls], **PLOT_CONFIG["bar_style"])
            bottom += frac

    # Direction labels
    ax.set_xticks(theta_positions)
    ax.set_xticklabels(DIRECTION_LABELS,
                       fontsize=PLOT_CONFIG["direction_label_fontsize"],
                       fontweight=PLOT_CONFIG["direction_label_fontweight"])

    # Radial scaling
    rmax = (max(sum(freq_data[d]) for d in DIRECTION_LABELS)
            if normalize == "global" else 1.0)
    rings = [p / 100.0 for p in PLOT_CONFIG["radial_percentages"]]
    ax.set_rmax(max(rmax, max(rings)) * 1.02)

    # Hide default y labels and add custom % labels
    ax.set_yticklabels([])
    label_angle = np.radians(PLOT_CONFIG["percent_label_angle_deg"])
    for pct in PLOT_CONFIG["radial_percentages"]:
        r = pct / 100.0
        if r <= ax.get_rmax():
            ax.text(label_angle, r, f"{pct}%",
                    fontsize=PLOT_CONFIG["percentage_label_fontsize"],
                    fontweight="bold", ha="center", va="bottom",
                    backgroundcolor="white", zorder=6)

    # Title
    param_full = "Zero-Crossing Period (Tz)" if parameter_name == "Tz" \
                 else "Peak Period (Tp)"
    title_norm = "Global frequency" if normalize == "global" else "Per-direction frequency"
    plt.title(f"Directional Distribution of Wave {param_full}\n{title_norm} by period class",
              fontsize=PLOT_CONFIG["title_fontsize"], fontweight="bold", pad=18)

    # Save
    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Windrose plot saved: {output_filename}")

def create_wave_period_legend(output_filename: str, fig_format: str, parameter_name: str) -> None:
    fig = plt.figure(figsize=LEGEND_CONFIG["figsize"])
    ax = fig.add_subplot(111)

    legend_elements = [
        Patch(facecolor=WAVE_PERIOD_COLORS[i], edgecolor="black", linewidth=LEGEND_CONFIG["edge_lw"])
        for i in range(len(WAVE_PERIOD_CLASS))
    ]

    param_full = "Tz" if parameter_name == "Tz" \
                 else "Peak Period (Tp)"
    legend_title = f"{param_full} Classes"

    font_props = font_manager.FontProperties(**LEGEND_CONFIG["font_properties"])
    leg = ax.legend(
        legend_elements, WAVE_PERIOD_CLASS,
        title=legend_title, loc="center", prop=font_props,
        frameon=True, fancybox=True, shadow=False, ncol=1,
    )
    leg.get_title().set_fontweight("bold")
    leg.get_title().set_fontsize(LEGEND_CONFIG["font_properties"]["size"] + 2)
    ax.axis("off")

    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Legend saved: {output_filename}")

# ==================== MAIN ====================

def main():
    print("=== WINDROSE PLOT GENERATOR FOR WAVE PERIOD (Tz/Tp) ===")
    try:
        args = parse_arguments()

        # Default formats now prefer JPG if unspecified (matches default filenames)
        plot_format  = _figure_format(args.output, args.format, default="jpg")
        legend_format = _figure_format(args.legend, args.legend_format, default="jpg")

        # Show working dir and absolute paths so you always know where files go
        print(f"[cwd] {os.path.abspath(os.getcwd())}")
        print(f"Input: {os.path.abspath(args.input)}")
        print(f"Output plot: {os.path.abspath(args.output)} (format={plot_format})")
        print(f"Output legend: {os.path.abspath(args.legend)} (format={legend_format})")
        print(f"Parameter: {args.parameter} | Normalize: {args.normalize}")

        validate_input_file(args.input)
        df = load_wave_data(args.input, args.sheet, args.parameter)

        directions = df["dir_w"].to_numpy(dtype=float)
        periods = df[args.parameter].to_numpy(dtype=float)

        # Basic stats (helpful log)
        print(f"Rows: {len(df)} | dir_w in [min={np.nanmin(directions):.1f}, max={np.nanmax(directions):.1f}]° | "
              f"{args.parameter} in [min={np.nanmin(periods):.2f}, max={np.nanmax(periods):.2f}] s")

        processed = process_wave_period_data(directions, periods, normalize=args.normalize)

        create_windrose_plot(processed, args.output, plot_format, args.parameter, args.normalize)
        create_wave_period_legend(args.legend, legend_format, args.parameter)

        print("\n=== COMPLETED SUCCESSFULLY ===")
        print(f"Directional sectors : {len(DIRECTION_LABELS)}")
        print(f"Period classes      : {len(WAVE_PERIOD_CLASS)}")

    except Exception as error:
        print(f"ERROR: {str(error)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
