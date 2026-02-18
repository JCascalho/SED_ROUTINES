"""
======================================================
Windrose plot generator for oceanographic current data
File: windrose_currents.py
======================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: November 2025

Description
-----------
This script generates a windrose-style plot representing the distribution of
near-bed ocean current speeds as a function of direction. Directions are binned
into 32 canonical compass sectors (11.25° resolution). Output figures are
publication-grade and can be exported in vector (PDF, SVG) or raster (PNG, JPG)
formats at high resolution.

Directional frequencies can be normalized either globally (relative to all
observations) or per direction, such that the frequency classes within each
sector sum to 100%.

Input
-----
- Excel file (default):
    INPUT_using_OUTPUT_CURRENTS_WAVES_DATA_BASE_WITH_DELTA_MIX.xlsx
- Sheet name (command-line argument --sheet / -s, default: "C0_C7")
- Required columns (can be overridden via command-line arguments):
    - Direction column: dir_c  (argument --direction-col)
    - Speed column:     u_z    (argument --speed-col)

Output
------
- Windrose plot file:
    (argument --output / -o, default: WIND_ROSE_dir_c_C0_C7.jpg)
- Legend file:
    (argument --legend / -l, default: legend_Currents.jpg)
- Output formats controlled by:
    --format, --legend-format  (pdf, png, svg, jpg)

Usage examples
--------------
Basic usage (output format inferred from filename):
  python S_7_4 - wind rose plot generator for current data.py \
      --sheet "C0_C7"

Explicit JPG output with absolute paths:
  python S_7_4 - wind rose plot generator for current data.py \
      -s "C0_C7" \
      -o "C:/Temp/WIND_ROSE_dir_c_C0_C7.jpg" --format jpg \
      --legend "C:/Temp/legend_Currents.jpg" --legend-format jpg

Per-direction normalization (each sector sums to 100%):
  python S_7_4 - wind rose plot generator for current data.py \
      -s "C0_C7" \
      -o "WIND_ROSE_dir_c_C0_C7.jpg" --format jpg \
      --normalize per-direction

The windrose visualization provides a compact and intuitive representation of
directional current forcing relevant to sediment transport interpretation.
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

# Current speed bins (m/s) and labels 
CURRENT_SPEED_BINS = [0, 0.10, 0.20, 0.30, 0.40, np.inf]
CURRENT_SPEED_CLASS_LABELS = [
    "<0.10 m/s",
    "0.10–0.20 m/s",
    "0.20–0.30 m/s",
    "0.30–0.40 m/s",
    "≥0.40 m/s",
]
CURRENT_SPEED_COLORS = [
    (0/255, 166/255, 255/255),   # <0.10 - blue
    (34/255, 127/255, 34/255),   # 0.10-0.20 - green
    (243/255, 231/255, 0/255),   # 0.20-0.30 - yellow
    (217/255, 154/255, 37/255),  # 0.30-0.40 orange
    (255/255, 0/255, 0/255),     # >0.40 - red
]

DIRECTION_BINS = np.linspace(0, 360, 33)
DIRECTION_LABELS = [
    "N", "NbE", "NNE", "NEbN", "NE", "NEbE", "ENE", "EbN",
    "E", "EbS", "ESE", "SEbE", "SE", "SEbS", "SSE", "SbE",
    "S", "SbW", "SSW", "SWbS", "SW", "SWbW", "WSW", "WbS",
    "W", "WbN", "WNW", "NWbW", "NW", "NWbN", "NNW", "NbW",
]

PLOT_CONFIG = {
    "figsize": (12, 12),
    "dpi": 600,
    "grid_style": {"linestyle": "-", "color": "k", "linewidth": 0.8, "alpha": 0.35},
    "bar_style": {"edgecolor": "black", "linewidth": 0.8, "alpha": 1.0},
    "radial_percentages": [2, 4, 6, 8],
    "direction_label_fontsize": 16,
    "direction_label_fontweight": "bold",
    "percentage_fontsize": 14,
    "percentage_fontweight": "bold",
    "percent_label_angle_deg": 120,
    "title_fontsize": 16,
}

LEGEND_CONFIG = {
    "figsize": (6, 3.0),
    "font_properties": {"weight": "bold", "size": 12},
    "title": "Current Speed Classes (m/s)",
    "edge_lw": 0.8,
}

# ==================== CLI ====================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality windrose plot for ocean current data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--sheet", "-s", default="C0_C7", #select the desired Excel file sheet 
                        help="Excel sheet name containing current data (default: C0_C7)") #select the desired Excel file sheet 
    parser.add_argument("--output", "-o", default="WIND_ROSE_dir_c_C0_C7.jpg", #select the desired Excel file sheet 
                        help="Output filename for windrose plot (default: WIND_ROSE_dir_c_C0_C7.jpg)") #select the desired Excel file sheet  
    parser.add_argument("--legend", "-l", default="legend_Currents.jpg",
                        help="Legend output filename (default: legend_Currents.jpg)")
    parser.add_argument("--input", "-i", default=INPUT_FILE,
                        help=f"Input Excel file path (default: {INPUT_FILE})")
    parser.add_argument("--speed-col", default="u_z",
                        help="Current speed column name (default: u_z)")
    parser.add_argument("--direction-col", default="dir_c",
                        help="Current direction column name (default: dir_c)")
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

def _sanitize_basename(name: str) -> str:
    safe = "".join((c if c not in _INVALID_CHARS else "_") for c in name)
    return safe.strip().strip(".")

def _trim_raw_format_suffix(name: str) -> str:
    """
    If the basename ends with a bare format word like 'jpg' (no preceding dot),
    strip it so we don't end up with '...jpg.jpg'.
    """
    for suf in ["jpg", "jpeg", "png", "pdf", "svg", "tif", "tiff", "webp"]:
        low = name.lower()
        if low.endswith(suf) and not low.endswith("." + suf):
            return name[: -len(suf)]
    return name

def _with_ext(basename: str, fig_format: str) -> str:
    fmt = (fig_format or "").lower()
    fmt = "jpg" if fmt in {"jpg", "jpeg"} else fmt
    root, ext = os.path.splitext(basename)
    wanted = "." + (fmt or "jpg")
    if ext.lower() != wanted:
        return root + wanted
    return basename

def _sanitize_and_fix_extension(path: str, fig_format: str) -> str:
    directory, base = os.path.split(path)
    # map ≥ ≤ > < first, then sanitize filename characters
    base = _semantic_token(base)
    base = _sanitize_basename(base)
    base = _trim_raw_format_suffix(base)
    base = _with_ext(base, fig_format)
    return os.path.join(directory if directory else "", base)

def _post_save_report(path: str) -> None:
    apath = os.path.abspath(path)
    print(f"[savefig] attempted: {apath}")
    print(f"[savefig] exists : {os.path.exists(apath)}")

def _print_supported_filetypes():
    ft = plt.gcf().canvas.get_supported_filetypes()
    print(f"[matplotlib] supported output formats: {', '.join(sorted(ft.keys()))}")

def _safe_savefig(path: str, fig_format: str, dpi: int, bbox_inches: str = "tight") -> None:
    fmt = (fig_format or "").lower()
    fmt = "jpeg" if fmt in {"jpg", "jpeg"} else fmt
    _ensure_parent_dir(path)
    try:
        plt.savefig(path, dpi=dpi, bbox_inches=bbox_inches, format=fmt if fmt else None)
    except ValueError:
        if fmt == "jpeg":
            fallback_path = os.path.splitext(path)[0] + ".png"
            print(f"Warning: JPEG export requires Pillow. Falling back to PNG: {fallback_path}", file=sys.stderr)
            plt.savefig(fallback_path, dpi=dpi, bbox_inches=bbox_inches, format="png")
            _post_save_report(fallback_path)
            return
        raise
    _post_save_report(path)

def _figure_format(filename: str, arg_format: str, default: str = "jpg") -> str:
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

def load_current_data(file_path: str, sheet_name: str, speed_column: str, direction_column: str) -> pd.DataFrame:
    xl = pd.ExcelFile(file_path)
    if sheet_name not in xl.sheet_names:
        raise ValueError(f"Sheet '{sheet_name}' not found. Available sheets: {xl.sheet_names}")
    df = pd.read_excel(file_path, sheet_name=sheet_name)

    required = [direction_column, speed_column]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Available: {list(df.columns)}")

    if (df[direction_column] < 0).any() or (df[direction_column] > 360).any():
        print("Warning: Some current directions outside 0–360°.", file=sys.stderr)
    if (df[speed_column] < 0).any():
        raise ValueError("Negative current speeds found.")

    high = df[speed_column] > 2.0
    if high.any():
        print(f"Warning: {int(high.sum())} current speeds > 2.0 m/s detected.", file=sys.stderr)

    start = len(df)
    df = df.dropna(subset=required)
    if df.empty:
        raise ValueError("No valid data remain after cleaning.")
    dropped = start - len(df)
    if dropped > 0:
        print(f"Note: dropped {dropped} rows with missing {required}.", file=sys.stderr)

    return df

# ==================== CORE LOGIC ====================

def classify_current_data(
    current_directions: np.ndarray,
    current_speeds: np.ndarray,
    normalize: str = "global"
) -> Dict[str, List[float]]:
    dir_idx = np.digitize(current_directions % 360.0, DIRECTION_BINS, right=False) - 1
    dir_idx = dir_idx % 32

    spd_idx = np.digitize(current_speeds, CURRENT_SPEED_BINS, right=False) - 1
    valid = (spd_idx >= 0) & (spd_idx < len(CURRENT_SPEED_BINS) - 1)

    counts = {lab: [0] * (len(CURRENT_SPEED_BINS) - 1) for lab in DIRECTION_LABELS}
    for d, s in zip(dir_idx[valid], spd_idx[valid]):
        counts[DIRECTION_LABELS[d]][s] += 1

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

def create_current_windrose_plot(
    freq_data: Dict[str, List[float]],
    output_filename: str,
    fig_format: str,
    normalize: str
) -> None:
    fig = plt.figure(figsize=PLOT_CONFIG["figsize"])
    ax = fig.add_subplot(111, polar=True)

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    theta_positions = np.linspace(0.0, 2 * np.pi, len(DIRECTION_LABELS), endpoint=False)
    sector_width = 2 * np.pi / len(DIRECTION_LABELS)

    ax.grid(True, **PLOT_CONFIG["grid_style"])

    for i, direction in enumerate(DIRECTION_LABELS):
        bottom = 0.0
        for cls in range(len(CURRENT_SPEED_CLASS_LABELS)):
            frac = freq_data[direction][cls]
            if frac <= 0:
                continue
            ax.bar(theta_positions[i], frac, width=sector_width, bottom=bottom,
                   color=CURRENT_SPEED_COLORS[cls], **PLOT_CONFIG["bar_style"])
            bottom += frac

    ax.set_xticks(theta_positions)
    ax.set_xticklabels(
        DIRECTION_LABELS,
        fontsize=PLOT_CONFIG["direction_label_fontsize"],
        fontweight=PLOT_CONFIG["direction_label_fontweight"],
    )

    rmax = (max(sum(freq_data[d]) for d in DIRECTION_LABELS)
            if normalize == "global" else 1.0)
    rings = [p / 100.0 for p in PLOT_CONFIG["radial_percentages"]]
    ax.set_rmax(max(rmax, max(rings)) * 1.02)

    ax.set_yticklabels([])
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

    title_norm = "Global frequency" if normalize == "global" else "Per-direction frequency"
    plt.title(f"Directional Distribution of Current Speed\n{title_norm} by speed class",
              fontsize=PLOT_CONFIG["title_fontsize"], fontweight="bold", pad=16)

    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Current windrose plot saved: {output_filename}")

def create_current_speed_legend(output_filename: str, fig_format: str) -> None:
    fig = plt.figure(figsize=LEGEND_CONFIG["figsize"])
    ax = fig.add_subplot(111)

    handles = [
        Patch(facecolor=CURRENT_SPEED_COLORS[i], edgecolor="black", linewidth=LEGEND_CONFIG["edge_lw"])
        for i in range(len(CURRENT_SPEED_CLASS_LABELS))
    ]

    font_props = font_manager.FontProperties(**LEGEND_CONFIG["font_properties"])
    leg = ax.legend(
        handles, CURRENT_SPEED_CLASS_LABELS,
        title=LEGEND_CONFIG["title"], loc="center", prop=font_props,
        frameon=True, fancybox=True, shadow=False, ncol=1,
    )
    leg.get_title().set_fontweight("bold")
    ax.axis("off")

    _print_supported_filetypes()
    _safe_savefig(output_filename, fig_format, dpi=PLOT_CONFIG["dpi"])
    plt.close()
    print(f"Current speed legend saved: {output_filename}")

# ==================== MAIN ====================

def main():
    print("=== WINDROSE PLOT GENERATOR FOR OCEANOGRAPHIC CURRENT DATA ===")
    try:
        args = parse_arguments()

        plot_format  = _figure_format(args.output, args.format, default="jpg")
        legend_format = _figure_format(args.legend, args.legend_format, default="jpg")

        # Apply semantic token mapping + sanitization + extension fix
        raw_output, raw_legend = args.output, args.legend
        args.output = _sanitize_and_fix_extension(args.output, plot_format)
        args.legend = _sanitize_and_fix_extension(args.legend, legend_format)

        print(f"[cwd] {os.path.abspath(os.getcwd())}")
        print(f"Input         : {os.path.abspath(args.input)}")
        if raw_output != args.output:
            print(f"Output plot   : {os.path.abspath(raw_output)}  →  {os.path.abspath(args.output)} (format={plot_format})")
        else:
            print(f"Output plot   : {os.path.abspath(args.output)} (format={plot_format})")
        if raw_legend != args.legend:
            print(f"Output legend : {os.path.abspath(raw_legend)}  →  {os.path.abspath(args.legend)} (format={legend_format})")
        else:
            print(f"Output legend : {os.path.abspath(args.legend)} (format={legend_format})")
        print(f"Speed column  : {args.speed_col} | Direction column: {args.direction_col} | Normalize: {args.normalize}")

        validate_input_file(args.input)
        df = load_current_data(args.input, args.sheet, args.speed_col, args.direction_col)

        directions = df[args.direction_col].to_numpy(dtype=float)
        speeds = df[args.speed_col].to_numpy(dtype=float)

        print(f"Rows: {len(df)} | {args.direction_col} in [min={np.nanmin(directions):.1f}, max={np.nanmax(directions):.1f}]° | "
              f"{args.speed_col} in [min={np.nanmin(speeds):.3f}, max={np.nanmax(speeds):.3f}] m/s")
        print(f"  > {100*(speeds >= 0.10).sum()/len(speeds):.1f}% ≥ 0.10 m/s (incipient motion)")
        print(f"  > {100*(speeds >= 0.20).sum()/len(speeds):.1f}% ≥ 0.20 m/s (significant transport)")

        processed = classify_current_data(directions, speeds, normalize=args.normalize)

        create_current_windrose_plot(processed, args.output, plot_format, args.normalize)
        create_current_speed_legend(args.legend, legend_format)

        print("\n=== COMPLETED SUCCESSFULLY ===")
        print(f"Directional sectors : {len(DIRECTION_LABELS)}")
        print(f"Speed classes       : {len(CURRENT_SPEED_BINS) - 1}")

        print("\nSEDIMENTOLOGICAL INTERPRETATION GUIDANCE:")
        print("- Dominant current directions indicate main transport pathways.")
        print("- Higher speed frequencies in specific sectors suggest energetic transport regimes.")
        print("- Seasonal variations can be analyzed by processing different time periods separately.")

    except Exception as error:
        print(f"ERROR: {str(error)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
