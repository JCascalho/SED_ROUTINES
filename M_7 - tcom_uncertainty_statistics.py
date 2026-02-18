"""
=======================================================================
Uncertainty statistics script for tracer centre-of-mass (TCOM) analysis
File:  tcom_uncertainty_statistics.py
======================================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: January 2026

This script implements an integrated workflow to quantify uncertainty in
tracer centre-of-mass (TCOM) positions derived from sediment tracer surveys.

ALL-IN-ONE functionality:
  1) Computes TCOM positional uncertainty using non-parametric bootstrap
     resampling in both Cartesian (x, y) and rotated along-shelf / cross-shelf
     (as, cs) coordinate systems.
       -> writes TRACER_CM_uncertainty_bootstrap.xlsx
  2) Generates combined multi-panel figures showing bootstrap uncertainty
     clouds and 95% covariance ellipses in relative axes (ΔE/ΔN).
       -> writes two PNG figures (C2–C4 and C5–C7)
  3) Exports ArcGIS-ready 95% confidence ellipse vertices to Excel
     (one worksheet per sampling campaign).
       -> writes CM_uncertainty_ellipses_95_vertices_by_campaign.xlsx

INPUT (expected in the same folder as this script):
  - input_TRACER_DATA.xlsx

-----------------------------------------------------------------------
IMPORTANT NOTE: Bootstrap realizations with zero tracer weight
-----------------------------------------------------------------------
Bootstrap resampling is performed by drawing sampling stations with
replacement and computing tracer-weighted centres of mass (CM) using
weights defined as TMcal = Ac / m.

For campaigns with patchy tracer recovery, some bootstrap realizations may
contain only stations with Ac = 0 (no detected tracer). In these cases,
the total tracer weight is zero (sum(TMcal) = 0) and the CM position is
undefined.

Such realizations do not contain physically meaningful information on tracer
CM position. Therefore, they are:
  - flagged as NaN during the bootstrap procedure, and
  - excluded prior to:
      (i) saving BOOT_C* sheets to Excel,
      (ii) computing percentile-based confidence intervals (CI2.5–CI97.5),
      (iii) estimating covariance matrices and confidence ellipses.

The number of valid and invalid bootstrap realizations is recorded for each
campaign in a dedicated Excel worksheet (BOOT_INVALID_COUNTS). This approach
preserves the integrity of the bootstrap procedure while ensuring that
summary statistics and uncertainty ellipses are computed exclusively from
physically meaningful TCOM realizations.
"""


from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS (edit here)
# ============================================================
INPUT_XLSX_NAME = "input_TRACER_DATA.xlsx"

# Outputs
BOOT_OUT_XLSX_NAME = "TRACER_CM_uncertainty_bootstrap.xlsx"
OUT_VERTICES_XLSX = "CM_uncertainty_ellipses_95_vertices_by_campaign.xlsx"

CAMPAIGNS_ALL = ["C2", "C3", "C4", "C5", "C6", "C7"]
N_PER_PAGE = 3

# Bootstrap settings
B = 20000
SEED = 42

# Station-position uncertainty (metres)
SIGMA_X_M = 2.5
SIGMA_Y_M = 2.5
SIGMA_AS_M = 2.5
SIGMA_CS_M = 2.5

# Figure settings (relative axes)
PAD_M = 5.0
ORIGIN_ROUND_M = 1.0   # set 10.0 if you want rounded origins
ELLIPSE_LW = 1.6

# ArcGIS vertex export density
NPTS_VERTICES = 180

# Original CM colors (your palette)
COLORS = {
    "C2": (255/255, 168/255,   0/255),  # orange
    "C3": (102/255,  51/255,   0/255),  # brown
    "C4": ( 22/255, 182/255,  78/255),  # green
    "C5": (  0/255, 231/255, 255/255),  # cyan
    "C6": (255/255,   0/255,   0/255),  # red
    "C7": (  0/255,  85/255, 255/255),  # blue
}
# ============================================================


# -------------------------
# Shared utilities
# -------------------------
def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"\n{label} not found:\n  {path}\n\n"
            f"Fix: Put '{path.name}' in the same folder as this script.\n"
        )

def get_ac_column(df: pd.DataFrame) -> str:
    """Robust detection of tracer mass column name."""
    if "Ac " in df.columns:
        return "Ac "
    if "Ac" in df.columns:
        return "Ac"
    raise KeyError("Neither 'Ac' nor 'Ac ' column found in the input workbook.")

def get_weights_TMcal(df: pd.DataFrame) -> np.ndarray:
    """Tracer weights: TMcal = Ac/m (safe for m<=0)."""
    ac_col = get_ac_column(df)
    if "m" not in df.columns:
        raise KeyError("Missing 'm' column.")
    ac = df[ac_col].to_numpy(float)
    m = df["m"].to_numpy(float)
    return np.where(m > 0, ac / m, 0.0)

def kish_n_eff(w: np.ndarray) -> float:
    """Kish effective sample size from normalized weights."""
    w = np.asarray(w, dtype=float)
    s = np.sum(w)
    if s <= 0 or np.isclose(s, 0.0):
        return np.nan
    p = w / s
    denom = np.sum(p**2)
    return float(1.0 / denom) if denom > 0 else np.nan


# -------------------------
# Part 1: Bootstrap computation + Excel output
# -------------------------
def weighted_cm(df: pd.DataFrame, w: np.ndarray) -> dict:
    """Deterministic CM for both coordinate systems."""
    for col in ["x", "y", "as", "cs"]:
        if col not in df.columns:
            raise KeyError(f"Missing '{col}' column.")

    s = float(np.sum(w))
    if s <= 0 or np.isclose(s, 0.0):
        return dict(CMx=np.nan, CMy=np.nan, CMas=np.nan, CMcs=np.nan, TM_sum=s)

    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)
    a = df["as"].to_numpy(float)
    c = df["cs"].to_numpy(float)

    return dict(
        CMx=float(np.sum(w * x) / s),
        CMy=float(np.sum(w * y) / s),
        CMas=float(np.sum(w * a) / s),
        CMcs=float(np.sum(w * c) / s),
        TM_sum=s,
    )

def bootstrap_cm_cloud(df: pd.DataFrame, B: int, seed: int,
                       sigma_x: float, sigma_y: float,
                       sigma_as: float, sigma_cs: float) -> pd.DataFrame:
    """
    Bootstrap CM positions by resampling stations with replacement and jittering coordinates.
    Any replicate with sum(weights)=0 is returned as NaN (invalid).
    """
    rng = np.random.default_rng(seed)

    w0 = get_weights_TMcal(df)
    x0 = df["x"].to_numpy(float)
    y0 = df["y"].to_numpy(float)
    as0 = df["as"].to_numpy(float)
    cs0 = df["cs"].to_numpy(float)

    n = len(df)
    out = np.full((B, 4), np.nan, dtype=float)

    for b in range(B):
        idx = rng.integers(0, n, size=n)
        w = w0[idx]
        s = np.sum(w)
        if s <= 0 or np.isclose(s, 0.0):
            # invalid replicate: CM undefined (all resampled stations had Ac=0)
            continue

        x = x0[idx] + rng.normal(0, sigma_x, size=n)
        y = y0[idx] + rng.normal(0, sigma_y, size=n)
        a = as0[idx] + rng.normal(0, sigma_as, size=n)
        c = cs0[idx] + rng.normal(0, sigma_cs, size=n)

        out[b, 0] = np.sum(w * x) / s
        out[b, 1] = np.sum(w * y) / s
        out[b, 2] = np.sum(w * a) / s
        out[b, 3] = np.sum(w * c) / s

    return pd.DataFrame(out, columns=["CMx", "CMy", "CMas", "CMcs"])

def summary_stats(v: np.ndarray) -> dict:
    v = np.asarray(v, dtype=float)
    return dict(
        Median=float(np.nanmedian(v)),
        Mean=float(np.nanmean(v)),
        CI2p5=float(np.nanpercentile(v, 2.5)),
        CI97p5=float(np.nanpercentile(v, 97.5)),
    )

def ellipse_params(x: np.ndarray, y: np.ndarray, level: float = 0.95) -> dict:
    """
    Covariance ellipse parameters:
      - a,b: semi-axes (metres)
      - angle_deg: degrees CCW from +x axis
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 3:
        return dict(a=np.nan, b=np.nan, angle_deg=np.nan)

    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    chi2 = 5.991 if abs(level - 0.95) < 1e-12 else 2.279
    a = float(np.sqrt(vals[0] * chi2))
    b = float(np.sqrt(vals[1] * chi2))
    angle = float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])))

    return dict(a=a, b=b, angle_deg=angle)

def compute_bootstrap_excel(here: Path) -> Path:
    """
    Runs the bootstrap analysis and writes BOOT_OUT_XLSX_NAME.
    Applies Option 1: drops invalid (NaN) bootstrap realizations before:
      - saving BOOT_C* sheets
      - computing summaries/ellipses
    Also writes BOOT_INVALID_COUNTS for transparency.
    """
    input_xlsx = here / INPUT_XLSX_NAME
    out_xlsx = here / BOOT_OUT_XLSX_NAME

    require_file(input_xlsx, "INPUT_XLSX")
    xls = pd.ExcelFile(input_xlsx)

    original_rows = []
    boot_median_rows = []
    summary_rows = []
    ellipse_xy_rows = []
    ellipse_ascs_rows = []
    boot_clouds = {}
    invalid_counts_rows = []

    rng = np.random.default_rng(SEED)

    for c in CAMPAIGNS_ALL:
        df = pd.read_excel(xls, sheet_name=c)
        w = get_weights_TMcal(df)
        neff = kish_n_eff(w)

        cm0 = weighted_cm(df, w)
        original_rows.append(dict(Campaign=c, n_eff=neff, **cm0))

        seed_c = int(rng.integers(0, 2**31 - 1))
        boot_df_raw = bootstrap_cm_cloud(
            df=df, B=B, seed=seed_c,
            sigma_x=SIGMA_X_M, sigma_y=SIGMA_Y_M,
            sigma_as=SIGMA_AS_M, sigma_cs=SIGMA_CS_M
        )

        # -------- OPTION 1: drop invalid realizations (all NaN rows) --------
        n_total = len(boot_df_raw)
        n_invalid = int(boot_df_raw.isna().any(axis=1).sum())
        boot_df = boot_df_raw.dropna(how="any").reset_index(drop=True)
        n_valid = len(boot_df)

        invalid_counts_rows.append({
            "Campaign": c,
            "B_requested": B,
            "B_total_rows": n_total,
            "B_invalid": n_invalid,
            "B_valid": n_valid,
            "p_invalid": (n_invalid / n_total) if n_total > 0 else np.nan
        })

        boot_clouds[c] = boot_df

        # Bootstrap median + CI (wide)
        row = {"Campaign": c, "n_eff": neff, "B_valid": n_valid, "B_invalid": n_invalid}
        for metric in ["CMx", "CMy", "CMas", "CMcs"]:
            stats = summary_stats(boot_df[metric].to_numpy(float))
            row[f"{metric}_median"] = stats["Median"]
            row[f"{metric}_CI2p5"] = stats["CI2p5"]
            row[f"{metric}_CI97p5"] = stats["CI97p5"]
        boot_median_rows.append(row)

        # Long summary
        for metric in ["CMx", "CMy", "CMas", "CMcs"]:
            stats = summary_stats(boot_df[metric].to_numpy(float))
            summary_rows.append(dict(Campaign=c, Metric=metric, n_eff=neff, B_valid=n_valid, **stats))

        # Ellipses (computed from VALID replicates only)
        ep_xy = ellipse_params(boot_df["CMx"], boot_df["CMy"], level=0.95)
        ellipse_xy_rows.append(dict(Campaign=c, level=0.95, B_valid=n_valid, **ep_xy))

        ep_as = ellipse_params(boot_df["CMas"], boot_df["CMcs"], level=0.95)
        ellipse_ascs_rows.append(dict(Campaign=c, level=0.95, B_valid=n_valid, **ep_as))

    df_original = pd.DataFrame(original_rows)
    df_boot_median = pd.DataFrame(boot_median_rows)
    df_summary = pd.DataFrame(summary_rows)
    df_ell_xy = pd.DataFrame(ellipse_xy_rows)
    df_ell_as = pd.DataFrame(ellipse_ascs_rows)
    df_invalid = pd.DataFrame(invalid_counts_rows)

    df_params = pd.DataFrame([dict(
        INPUT_XLSX=str(input_xlsx),
        B_requested=B,
        SEED=SEED,
        SIGMA_X_M=SIGMA_X_M,
        SIGMA_Y_M=SIGMA_Y_M,
        SIGMA_AS_M=SIGMA_AS_M,
        SIGMA_CS_M=SIGMA_CS_M,
        NOTE="Weights TMcal = Ac/m; invalid bootstrap replicates (sum weights=0) removed before summaries/ellipses."
    )])

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        df_params.to_excel(writer, sheet_name="PARAMS", index=False)
        df_invalid.to_excel(writer, sheet_name="BOOT_INVALID_COUNTS", index=False)

        df_original.to_excel(writer, sheet_name="ORIGINAL_CM", index=False)
        df_boot_median.to_excel(writer, sheet_name="BOOT_MEDIAN_CM", index=False)
        df_summary.to_excel(writer, sheet_name="COM_SUMMARY", index=False)
        df_ell_xy.to_excel(writer, sheet_name="ELLIPSE95_XY", index=False)
        df_ell_as.to_excel(writer, sheet_name="ELLIPSE95_ASCS", index=False)

        for c in CAMPAIGNS_ALL:
            boot_clouds[c].to_excel(writer, sheet_name=f"BOOT_{c}"[:31], index=False)

    print("Saved bootstrap Excel:", out_xlsx.name)
    return out_xlsx


# -------------------------
# Part 2: Postprocess figures + ArcGIS vertices
# -------------------------
def ellipse_points(cx, cy, a, b, angle_deg, n=240):
    """Polyline points for plotting."""
    if not np.isfinite(a) or not np.isfinite(b) or a <= 0 or b <= 0:
        return None
    t = np.linspace(0, 2*np.pi, n)
    ex = a * np.cos(t)
    ey = b * np.sin(t)
    ang = np.deg2rad(angle_deg)
    R = np.array([[np.cos(ang), -np.sin(ang)],
                  [np.sin(ang),  np.cos(ang)]])
    xy = R @ np.vstack((ex, ey))
    x = xy[0] + cx
    y = xy[1] + cy
    return x, y

def ellipse_vertices(cx, cy, a, b, angle_deg, npts=180):
    """Vertex list for ArcGIS polygon construction."""
    t = np.linspace(0, 2*np.pi, npts, endpoint=False)
    x = a * np.cos(t)
    y = b * np.sin(t)
    ang = np.deg2rad(angle_deg)
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return np.column_stack([cx + xr, cy + yr])

def compute_page_origin_and_limits(dfs: list[pd.DataFrame]) -> tuple[float, float, float, float, float, float]:
    """Relative-axes origin + shared square window."""
    all_x = np.concatenate([df["x"].to_numpy(float) for df in dfs])
    all_y = np.concatenate([df["y"].to_numpy(float) for df in dfs])
    xmin, xmax = np.nanmin(all_x) - PAD_M, np.nanmax(all_x) + PAD_M
    ymin, ymax = np.nanmin(all_y) - PAD_M, np.nanmax(all_y) + PAD_M
    side = max(xmax - xmin, ymax - ymin)
    xmid, ymid = 0.5*(xmin + xmax), 0.5*(ymin + ymax)
    xmin, xmax = xmid - side/2, xmid + side/2
    ymin, ymax = ymid - side/2, ymid + side/2
    if ORIGIN_ROUND_M > 0:
        x0 = np.floor(xmin / ORIGIN_ROUND_M) * ORIGIN_ROUND_M
        y0 = np.floor(ymin / ORIGIN_ROUND_M) * ORIGIN_ROUND_M
    else:
        x0, y0 = xmin, ymin
    return x0, y0, xmin, xmax, ymin, ymax

def make_page(fig_path: Path, xls_in: pd.ExcelFile,
              df_original: pd.DataFrame, df_boot_median: pd.DataFrame, df_ell_xy: pd.DataFrame,
              campaigns: list[str]):

    dfs = [pd.read_excel(xls_in, sheet_name=c) for c in campaigns]
    x0, y0, xmin, xmax, ymin, ymax = compute_page_origin_and_limits(dfs)

    ncols = len(campaigns)
    fig, axes = plt.subplots(3, ncols, figsize=(4.6 * ncols, 10),
                             gridspec_kw=dict(height_ratios=[1, 1, 0.6]))

    for i, (c, df) in enumerate(zip(campaigns, dfs)):
        col = COLORS.get(c, (0.2, 0.2, 0.2))

        ac_col = get_ac_column(df)
        ac = df[ac_col].to_numpy(float)
        has_tracer = ac > 0

        orig_row = df_original[df_original["Campaign"] == c].iloc[0]
        med_row = df_boot_median[df_boot_median["Campaign"] == c].iloc[0]
        ell_row = df_ell_xy[df_ell_xy["Campaign"] == c].iloc[0]

        CMx, CMy = float(orig_row["CMx"]), float(orig_row["CMy"])
        n_eff = float(orig_row["n_eff"])

        new_x, new_y = float(med_row["CMx_median"]), float(med_row["CMy_median"])
        a, b, ang = float(ell_row["a"]), float(ell_row["b"]), float(ell_row["angle_deg"])

        # Row 1: geometry
        ax = axes[0, i]
        ax.scatter(df.loc[~has_tracer, "x"] - x0, df.loc[~has_tracer, "y"] - y0,
                   facecolors="none", edgecolors="0.6", s=26, linewidths=0.8, zorder=1)
        ax.scatter(df.loc[has_tracer, "x"] - x0, df.loc[has_tracer, "y"] - y0,
                   c="black", s=30, zorder=2)
        ax.scatter(CMx - x0, CMy - y0, c=[col], s=75, edgecolors="black", linewidths=0.7, zorder=3)
        ax.set_title(c)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(xmin - x0, xmax - x0)
        ax.set_ylim(ymin - y0, ymax - y0)

        # Row 2: uncertainty
        ax = axes[1, i]
        ax.scatter(df["x"] - x0, df["y"] - y0,
                   facecolors="white", edgecolors="black", s=18, linewidths=0.55, zorder=1)
        pts = ellipse_points(new_x - x0, new_y - y0, a, b, ang, n=240)
        if pts is not None:
            ax.plot(pts[0], pts[1], color="black", lw=ELLIPSE_LW, zorder=2)
        ax.scatter(new_x - x0, new_y - y0, marker="s", c="0.55", s=75,
                   edgecolors="black", linewidths=0.7, zorder=3)
        ax.scatter(CMx - x0, CMy - y0, c=[col], s=65,
                   edgecolors="black", linewidths=0.7, zorder=4)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(xmin - x0, xmax - x0)
        ax.set_ylim(ymin - y0, ymax - y0)

        # Row 3: n_eff
        ax = axes[2, i]
        ax.bar([0], [n_eff], color=col, width=0.6)
        ax.set_ylim(0, len(df))
        ax.set_xticks([])
        if i == 0:
            ax.set_ylabel("n_eff")

    axes[0, 0].set_ylabel("Sampling geometry")
    axes[1, 0].set_ylabel("Tracer uncertainty")
    fig.supxlabel(f"ΔEasting (m) from {int(round(x0))} m — ETRS89 / UTM 29N (EPSG:25829)", y=0.03)
    fig.supylabel(f"ΔNorthing (m) from {int(round(y0))} m — ETRS89 / UTM 29N (EPSG:25829)", x=0.02)
    fig.suptitle("Sampling geometry, tracer-weighted uncertainty (95% ellipse), and effective sample size",
                 fontsize=14, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)

def make_all_pages(here: Path, xls_in: pd.ExcelFile,
                   df_original: pd.DataFrame, df_boot_median: pd.DataFrame, df_ell_xy: pd.DataFrame):
    pages = [CAMPAIGNS_ALL[i:i + N_PER_PAGE] for i in range(0, len(CAMPAIGNS_ALL), N_PER_PAGE)]
    for camps in pages:
        tag = f"{camps[0]}-{camps[-1]}"
        out = here / f"Combined_geometry_uncertainty_neff_{tag}_relativeAxes.png"
        make_page(out, xls_in, df_original, df_boot_median, df_ell_xy, camps)
        print("Saved figure:", out.name)

def export_vertices_by_campaign(here: Path, boot_xlsx: Path):
    df_cm = pd.read_excel(boot_xlsx, sheet_name="BOOT_MEDIAN_CM")
    df_ell = pd.read_excel(boot_xlsx, sheet_name="ELLIPSE95_XY")
    df = df_cm.merge(df_ell, on="Campaign", how="inner").set_index("Campaign")

    out_path = here / OUT_VERTICES_XLSX
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for c in CAMPAIGNS_ALL:
            if c not in df.index:
                continue
            r = df.loc[c]
            verts = ellipse_vertices(
                r["CMx_median"], r["CMy_median"],
                r["a"], r["b"], r["angle_deg"],
                npts=NPTS_VERTICES
            )
            out_df = pd.DataFrame({
                "Order": np.arange(len(verts), dtype=int),
                "X": verts[:, 0].astype(float),
                "Y": verts[:, 1].astype(float),
            })
            out_df.to_excel(writer, sheet_name=str(c)[:31], index=False)

        meta = pd.DataFrame([{
            "source_bootstrap_excel": str(boot_xlsx),
            "npts_per_ellipse": NPTS_VERTICES,
            "crs": "ETRS89 / UTM Zone 29N (EPSG:25829)",
            "notes": "Ellipse vertices computed from (CMx_median, CMy_median, a, b, angle_deg)."
        }])
        meta.to_excel(writer, sheet_name="META", index=False)

    print("Saved ArcGIS vertices Excel:", out_path.name)


def main():
    here = Path(__file__).resolve().parent

    # 1) Compute bootstrap Excel (always)
    boot_xlsx = compute_bootstrap_excel(here)

    # 2) Postprocess using freshly created Excel
    input_xlsx = here / INPUT_XLSX_NAME
    xls_in = pd.ExcelFile(input_xlsx)

    df_original = pd.read_excel(boot_xlsx, sheet_name="ORIGINAL_CM")
    df_boot_median = pd.read_excel(boot_xlsx, sheet_name="BOOT_MEDIAN_CM")
    df_ell_xy = pd.read_excel(boot_xlsx, sheet_name="ELLIPSE95_XY")

    make_all_pages(here, xls_in, df_original, df_boot_median, df_ell_xy)
    export_vertices_by_campaign(here, boot_xlsx)

    print("\nDONE. Outputs written next to this script.")


if __name__ == "__main__":
    main()
