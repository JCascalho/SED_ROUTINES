"""
============================================
Windrose plot generator for wave height data
File: windrose_wave_height.py
============================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Description
-----------
This script generates a windrose-style plot showing the directional
distribution of significant wave height (Hs). Directions are binned into
32 canonical compass sectors (11.25° resolution). Output figures are
publication-grade and can be exported in vector (PDF, SVG) or raster
(PNG, JPG) formats at high resolution.

The resulting visualization provides a concise overview of directional
wave forcing relevant to sediment mobility and tracer transport analyses.

Input
-----
- Excel file (default):
    INPUT_using_OUTPUT_CURRENTS_WAVES_DATA_BASE_WITH_DELTA_MIX.xlsx
- Sheet name (command-line argument --sheet / -s, default: "C0_C7")
- Required columns:
    - dir_w : wave direction in degrees (0–360; “from” convention)
    - Hs    : significant wave height [m]

Output
------
- Windrose plot file:
    (argument --output / -o, default: WINDROSE_Hs_C0_C7.jpg)
- Legend file:
    (argument --legend / -l, default: legend_Hs.jpg)
- Output formats controlled by:
    --format, --legend-format  (pdf, png, svg, jpg)

Usage examples
--------------
Basic usage (output format inferred from filename):
  python S_7_5 - wind rose plot generator for wave height data.py \
      --sheet "C0_C7_tau_max>0.000"

Explicit JPG output with absolute paths:
  python S_7_5 - wind rose plot generator for wave height data.py \
      -s "C0_C7_tau_max>0.000" \
      -o "C:/Temp/WINDROSE_Hs.jpg" --format jpg \
      --legend "C:/Temp/legend_Hs.jpg" --legend-format jpg

Per-direction normalization (each sector sums to 100%):
  python S_7_5 - wind rose plot generator for wave height data.py \
      -s "C0_C7_tau_max>0.000" \
      -o "WINDROSE_Hs.jpg" --format jpg --normalize per-direction

Note (Windows shell)
-------------------
If the sheet name contains the '>' character (e.g., tau_max>0.000), it must
be quoted on the command line:
    --sheet "C0_C7_tau_max>0.000"

Frequencies can be normalized either globally (relative to all observations)
or per direction, such that the frequency classes within each sector sum to
100%.
"""
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

# Wave height bins (m) and colors
WAVE_HEIGHT_BINS = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, np.inf]
WAVE_HEIGHT_COLORS = [
    (0/255, 166/255, 255/255),   # very low
    (0/255, 183/255, 166/255),   # low
    (255/255, 242/255, 0/255),   # moderate
    (255/255, 161/255, 0/255),   # medium
    (255/255, 0/255, 0/255),     # high
    (153/255, 102/255, 51/255),  # very high
    (51/255, 51/255, 51/255),    # extreme
]
# Sanity check: colors must match number of classes
N_CLASSES = len(WAVE_HEIGHT_BINS) - 1
assert len(WAVE_HEIGHT_COLORS) == N_CLASSES, "Colors must match number of Hs classes"

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
    "radial_percentages": [5, 10, 15, 20],  # adapt the reference rings according to the maximum % value 
    "direction_label_fontsize": 16,
    "direction_label_fontweight": "bold",
    "percentage_fontsize": 14,
    "percentage_fontweight": "bold",
    "percent_label_angle_deg": 85,  # angle where % labels are placed
}

LEGEND_CONFIG = {
    "figsize": (6, 3.2),
    "font_properties": {"weight": "bold", "size": 12},
    "title": "Hs classes (m)",
    "edge_lw": 0.6,
}

# ==================== CLI ====================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality windrose plot for wave height (Hs).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --sheet "C0_C2_tau_max>0.000"
  %(prog)s -s "C3_C4_tau_max>0.259" -o results/wave_analysis.jpg --format jpg
""",
    )
    parser.add_argument("--sheet", "-s", default="C0_C7", #select the desired Excel file sheet 
                        help="Excel sheet name to process (default: C0_C7)") #select the desired Excel file sheet
    # Default to JPG filenames so Spyder/no-args users get JPGs by default
    parser.add_argument("--output", "-o", default="WIND_ROSE_Hs_C0_C7.jpg", 
                        help="Output filename for windrose plot (default: WIND_ROSE_Hs_C0_C7.jpg)")
    parser.add_argument("--legend", "-l", default="legend_Hs.jpg",
                        help="Legend output filename (default: legend_Hs.jpg)")
    parser.add_argument("--input", "-i", default=INPUT_FILE,
                        help=f"Input Excel file path (default: {INPUT_FILE})")
    parser.add_argument("--format", choices=["pdf", "png", "svg", "jpg"], default=None,
                        help="Override plot format (inferred from filename if omitted)")
    parser.add_argument("--legend-format", choices=["pdf", "png", "svg", "jpg"], default=None,
                        help="Override legend format (inferred from filename if omitted)")
    parser.add_argument("--normalize", choices=["global", "per-direction"], default="global",
                        help="Normalization: 'global' (default) or 'per-direction'")
    return parser.parse_args()

# ==================== SAVEFIG HELPERS (semantic tokens + robust JPG) ====================

_INVALID_CHARS = '<>:"/\\|?*'

def _semantic_token(s: str) -> str:
    """Map relational symbols to readable tokens before sanitization."""
    return (
        s.replace("≥", "_ge_")
         .replace("≤", "_le_")
         .replace(">", "_gt_")
         .replace("<", "_lt_")
    )

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
        # Use the exception object so linters don't complain and you get context
        print(f"[savefig] ValueError during save: {e}", file=sys.stderr)
        raise
    _post_save_report(path)

def _figure_format(filename: str, arg_format: str, default: str = "jpg") -> str:
    # default changed to 'jpg' so even with no args we try to write JPG
    ext = os.path.splitext(filename)[1].lower().lstrip(".")
    if arg_format:
        return arg_format
    if ext in {"pdf", "png", "svg", "jpg", "jpeg"}:
        return "jpg" if ext == "jpeg" else ext
    return default

# ==================== I/O & VALIDATION ====================

def validate_file(file_path: str) -> None:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"No read permission for file: {file_path}")

def load_data(file_path: str, sheet_name: str) -> pd.DataFrame:
    xl = pd.ExcelFile(file_path)
    if sheet_name not in xl.sheet_names:
        raise ValueError(f"Sheet '{sheet_name}' not found. Available: {xl.sheet_names}")
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    required = ["dir_w", "Hs"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Available: {list(df.columns)}")

    # Basic checks
    if (df["dir_w"] < 0).any() or (df["dir_w"] > 360).any():
        print("Warning: Some 'dir_w' values outside 0–360°.", file=sys.stderr)
    if (df["Hs"] < 0).any():
        raise ValueError("Negative Hs values found.")

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

def process_wave_data(
    wave_directions: np.ndarray,
    wave_heights: np.ndarray,
    normalize: str = "global"
) -> Dict[str, List[float]]:
    """Bin directions and Hs, then normalize counts."""
    # Wrap to [0,360) and bin robustly; ensure 0° -> first sector
    dir_idx = np.digitize(wave_directions % 360.0, DIRECTION_BINS, right=False) - 1
    dir_idx = dir_idx % 32  # handle -1 and 32 safely

    hs_idx = np.digitize(wave_heights, WAVE_HEIGHT_BINS, right=False) - 1
    valid = (hs_idx >= 0) & (hs_idx < len(WAVE_HEIGHT_BINS) - 1)

    counts = {lab: [0] * (len(WAVE_HEIGHT_BINS) - 1) for lab in DIRECTION_LABELS}
    for d, h in zip(dir_idx[valid], hs_idx[valid]):
        counts[DIRECTION_LABELS[d]][h] += 1

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
    processed: Dict[str, List[float]],
    output_filename: str,
    fig_format: str,
    normalize: str
) -> None:
    fig = plt.figure(figsize=PLOT_CONFIG["figsize"])
    ax = fig.add_subplot(111, polar=True)

    # Oceanographic orientation: 0° at North, clockwise
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    theta = np.linspace(0.0, 2 * np.pi, len(DIRECTION_LABELS), endpoint=False)
    width = 2 * np.pi / len(DIRECTION_LABELS)

    # Grid
    ax.grid(True, which="both", axis="both", **PLOT_CONFIG["grid_style"])

    # Stacked bars
    for i, direction in enumerate(DIRECTION_LABELS):
        bottom = 0.0
        for j, frac in enumerate(processed[direction]):
            if frac <= 0:
                continue
            ax.bar(theta[i], frac, width=width, bottom=bottom,
                   color=WAVE_HEIGHT_COLORS[j], **PLOT_CONFIG["bar_style"])
            bottom += frac

    # Direction labels
    ax.set_xticks(theta)
    ax.set_xticklabels(
        DIRECTION_LABELS,
        fontsize=PLOT_CONFIG["direction_label_fontsize"],
        fontweight=PLOT_CONFIG["direction_label_fontweight"],
    )

    # Radial scaling
    rmax = (max(sum(processed[d]) for d in DIRECTION_LABELS)
            if normalize == "global" else 1.0)
    rings = [p / 100.0 for p in PLOT_CONFIG["radial_percentages"]]
    ax.set_rmax(max(rmax, max(rings)) * 1.02)

    # Custom percentage labels
    ax.set_yticklabels([])  # hide defaults
    label_angle = np.radians(PLOT_CONFIG["percent_label_angle_deg"])
    for pct in PLOT_CONFIG["radial_percentages"]:
        r = pct / 100.0
        if r <= ax.get_rmax():
            ax.text(
                label_angle, r, f"{pct}",
                fontsize=PLOT_CONFIG["percentage_fontsize"],
                fontweight=PLOT_CONFIG["percentage_fontweight"],
                ha="center", va="bottom", backgroundcolor="white", zorder=6
            )

    # Title
    title_norm = "Global frequency" if normalize == "global" else "Per-direction frequency"
    plt.title(f"Directional Distribution of Wave Height (Hs)\n{title_norm} by Hs class",
              fontsize=16, fontweight="bold", pad=16)

    # Show supported types (to confirm 'jpeg' is active) and save
    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Windrose plot saved to: {output_filename}")

def create_legend(output_filename: str, fig_format: str) -> None:
    fig = plt.figure(figsize=LEGEND_CONFIG["figsize"])
    ax = fig.add_subplot(111)

    handles = [
        Patch(facecolor=WAVE_HEIGHT_COLORS[i], edgecolor="k", linewidth=LEGEND_CONFIG["edge_lw"])
        for i in range(len(WAVE_HEIGHT_COLORS))
    ]

    # Legend labels derived from bins
    labels = [f"< {WAVE_HEIGHT_BINS[1]}"] + [
        (f"{WAVE_HEIGHT_BINS[i]} – {WAVE_HEIGHT_BINS[i+1]}" if i < len(WAVE_HEIGHT_BINS) - 2
         else f"≥ {WAVE_HEIGHT_BINS[i]}")
        for i in range(1, len(WAVE_HEIGHT_BINS) - 1)
    ]

    font_props = font_manager.FontProperties(**LEGEND_CONFIG["font_properties"])
    leg = ax.legend(
        handles, labels,
        title=LEGEND_CONFIG["title"],
        loc="center",
        prop=font_props,
        frameon=True,
        fancybox=True,
        shadow=False,
        ncol=1,
    )
    leg.get_title().set_fontweight("bold")
    ax.axis("off")

    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Legend saved to: {output_filename}")

# ==================== MAIN ====================

def main():
    print("=== WINDROSE PLOT GENERATOR FOR Hs ===")
    print("32-point compass labels (canonical), publication-ready output\n")
    try:
        args = parse_arguments()

        # Default formats now prefer JPG if unspecified (matches default filenames)
        plot_format  = _figure_format(args.output, args.format, default="jpg")
        legend_format = _figure_format(args.legend, args.legend_format, default="jpg")

        # Show working dir so you know where relative paths go
        print(f"[cwd] {os.path.abspath(os.getcwd())}")
        print(f"Input: {os.path.abspath(args.input)}")
        print(f"Output plot: {os.path.abspath(args.output)} (format={plot_format})")
        print(f"Output legend: {os.path.abspath(args.legend)} (format={legend_format})")

        validate_file(args.input)
        df = load_data(args.input, args.sheet)

        dirs = df["dir_w"].to_numpy(dtype=float)
        hs = df["Hs"].to_numpy(dtype=float)

        print(f"Rows: {len(df)} | dir_w in [min={np.nanmin(dirs):.1f}, max={np.nanmax(dirs):.1f}]° | "
              f"Hs in [min={np.nanmin(hs):.2f}, max={np.nanmax(hs):.2f}] m")

        processed = process_wave_data(dirs, hs, normalize=args.normalize)

        create_windrose_plot(processed, args.output, plot_format, args.normalize)
        create_legend(args.legend, legend_format)

        print("\n=== COMPLETED SUCCESSFULLY ===")
        print(f"Directional sectors : {len(DIRECTION_LABELS)}")
        print(f"Hs classes          : {len(WAVE_HEIGHT_BINS) - 1}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
