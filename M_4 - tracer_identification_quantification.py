"""
==================================================================================
Tracer identification and quantification script
File: tracer_identification_quantification.py
==================================================================================

Author: João Cascalho (FCUL / IDL – University of Lisbon)
Date: January 2026

Tracer identification and quantification are performed using the System of
Automatic N-particle Detection (SAND), an automated tracer-particle detection
approach originally developed by Silva (2005) in the MATLAB environment using
the Image Processing Toolbox. The method was subsequently applied by Silva
et al. (2007) and Bosnic et al. (2017), with minor adaptations for different
tracer types and imaging conditions.

The SAND workflow is based on color classification in the CIE Lab* color space.
Digital images are converted from RGB to Lab*, and fluorescent tracer pixels
are identified using fixed threshold values applied to the L*, a*, and b*
channels to distinguish tracer grains from the natural sediment background.
The threshold values implemented here reproduce those used in the original
MATLAB SedPhoto workflow.

The resulting binary tracer masks are processed using a sequence of
morphological operations (erosion and dilation) designed to suppress
background noise while preserving small and low-intensity fluorescent
particles. Individual tracer grains are then identified as connected
components within the processed masks. Tracer abundance is quantified both
as the number of detected tracer particles and as tracer concentration,
expressed as the proportion of tracer pixels relative to the total image area.

In the present study, the original MATLAB-based image-analysis workflow was
reimplemented in Python using equivalent color-space transformations,
thresholding criteria, and morphological filtering steps. Validation against
the MATLAB implementation showed near-identical tracer counts across all
samples (differences < 1%). The Python implementation additionally
demonstrated improved robustness under low tracer-density conditions,
reliably detecting isolated tracer grains that are difficult to identify
visually or using less stringent image-processing approaches.

In addition to tracer counts, the area of individual tracer grains is
computed from the binary masks using connected-component analysis. Particle
areas are calculated in pixel units and converted to physical units (mm²)
using a calibrated pixel-to-area relationship, enabling subsequent grain-
scale analyses.

This script implements both tracer-particle detection and grain-area
quantification. It produces:
  (i)  per-image tracer counts and concentrations,
  (ii) binary tracer masks and outlined diagnostic images, and
  (iii) tables of individual particle areas derived from the processed masks.
"""


import os
import cv2
import numpy as np
import pandas as pd


# ============================================================
# SETTINGS — PART A (Tracer count and analysis)
# ============================================================

# MATLAB thresholds (uint8 Lab space, 0..255)
LORANGE = 35
AORANGE = 80
BORANGE = 130

# Morphology (SedPhoto BWFilter)
DISK_RADIUS = 2

# Output folders/files (keep your current names if you want)
MASK_DIR = "mask"
OUTLINE_DIR = "outlined"
EXCEL_TRACERCOUNT = "results_tracer_counts_C6_sample_26.xlsx" # use the appropriate name

# Image extensions (input photos)
IMG_EXTS = (".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp")


# ============================================================
# SETTINGS — PART B (Particle areas from masks)
# ============================================================

EXCEL_AREAS = "particle_areas_all_images_mm2_C6_sample_26.xlsx" # use the appropriate name

# Calibration
PIXEL_AREA_MM2 = 4.1698e-03

# Robust binarization for masks
BIN_THRESH = 128          # critical (do NOT use >0)
CONNECTIVITY = 8
MIN_AREA_PX = 0           # set 2–10 if you ever see tiny speckles


# ============================================================
# Helpers — PART A
# ============================================================

def disk_kernel(radius: int) -> np.ndarray:
    """Binary disk kernel similar to MATLAB strel('disk', r)."""
    r = int(radius)
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    k = ((x * x + y * y) <= r * r).astype(np.uint8)
    return k


def bwfilter_sedphoto(bw_255: np.ndarray, se: np.ndarray) -> np.ndarray:
    """MATLAB SedPhoto BWFilter: erode -> dilate -> dilate -> erode."""
    bw = bw_255
    bw = cv2.erode(bw, se, iterations=1)
    bw = cv2.dilate(bw, se, iterations=1)
    bw = cv2.dilate(bw, se, iterations=1)
    bw = cv2.erode(bw, se, iterations=1)
    return bw


def segment_sedphoto(bgr: np.ndarray, se: np.ndarray) -> np.ndarray:
    """Segment tracer pixels in uint8 Lab and apply BWFilter."""
    lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2LAB)
    L, a, b = cv2.split(lab)
    bw = ((L > LORANGE) & (a > AORANGE) & (b > BORANGE)).astype(np.uint8) * 255
    bw = bwfilter_sedphoto(bw, se)
    return bw


def count_particles_bwboundaries_like(bw_255: np.ndarray) -> int:
    """Approx MATLAB bwboundaries(BW,'noholes'): count external contours only."""
    contours, _ = cv2.findContours(bw_255, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    return int(len(contours))


def outline_image(bgr: np.ndarray, bw_255: np.ndarray) -> np.ndarray:
    """Draw external boundaries in white."""
    out = bgr.copy()
    contours, _ = cv2.findContours(bw_255, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(out, contours, -1, (255, 255, 255), 1)
    return out


def mask_filename_png(image_filename: str) -> str:
    """
    Always save masks as PNG (lossless), keeping 'mask_' prefix.
    Example: IMG_7023.JPG -> mask_IMG_7023.png
    """
    base = os.path.splitext(os.path.basename(image_filename))[0]
    return f"mask_{base}.png"


# ============================================================
# Part A — Tracer count and analysis
# ============================================================

def run_tracer_count_and_analysis():
    print("Working directory:", os.getcwd())

    files = sorted([f for f in os.listdir(".") if f.lower().endswith(IMG_EXTS)])
    if not files:
        raise RuntimeError("No image files found in the current working directory.")

    os.makedirs(MASK_DIR, exist_ok=True)
    os.makedirs(OUTLINE_DIR, exist_ok=True)

    se = disk_kernel(DISK_RADIUS)

    rows = []
    for fn in files:
        print("Processing image:", fn)

        bgr = cv2.imread(fn, cv2.IMREAD_COLOR)
        if bgr is None:
            print("  ⚠️ skipped (cannot read)")
            continue

        h, w = bgr.shape[:2]
        pix_photo = int(h * w)

        bw = segment_sedphoto(bgr, se)

        pix_particles = int(np.sum(bw > 0))
        concentration = float(pix_particles / pix_photo)
        nparts = count_particles_bwboundaries_like(bw)

        # Save mask as PNG (lossless)
        mask_name = mask_filename_png(fn)
        cv2.imwrite(os.path.join(MASK_DIR, mask_name), bw)

        # Save outlined image (same extension as input)
        outlined = outline_image(bgr, bw)
        cv2.imwrite(os.path.join(OUTLINE_DIR, f"outlined_{fn}"), outlined)

        rows.append({
            "image_name": fn,
            "np": nparts,
            "concentration": concentration,
            "pix_particles": pix_particles,
            "pix_photo": pix_photo
        })

        print(f"  np={nparts}  concentration={concentration:.6f}")

    df = pd.DataFrame(rows)
    df.to_excel(EXCEL_TRACERCOUNT, index=False)

    print("\nTracer-count analysis saved:", EXCEL_TRACERCOUNT)
    print("Masks saved in:", MASK_DIR)
    print("Outlined images saved in:", OUTLINE_DIR)


# ============================================================
# Part B — Particle areas from masks
# ============================================================

def clean_image_name_from_mask(mask_fname: str) -> str:
    """
    Convert:
      mask_IMG_7023.png -> IMG_7023
    """
    name = mask_fname
    if name.lower().startswith("mask_"):
        name = name[5:]
    name = os.path.splitext(name)[0]
    return name


def run_particle_areas_from_masks():
    if not os.path.isdir(MASK_DIR):
        raise RuntimeError(f"Mask folder not found: {MASK_DIR}")

    mask_files = sorted([
        f for f in os.listdir(MASK_DIR)
        if f.lower().endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"))
    ])
    if not mask_files:
        raise RuntimeError("No mask files found in mask/")

    rows = []
    for fname in mask_files:
        path = os.path.join(MASK_DIR, fname)
        img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        if img is None:
            print("Skipping unreadable mask:", fname)
            continue

        # Robust binarization (prevents speckle artifacts)
        bw01 = (img >= BIN_THRESH).astype(np.uint8)

        n, labels, stats, _ = cv2.connectedComponentsWithStats(
            bw01, connectivity=CONNECTIVITY
        )

        # Renumber particle_id sequentially after filtering
        pid = 0
        for lab in range(1, n):
            area_px = int(stats[lab, cv2.CC_STAT_AREA])
            if area_px < MIN_AREA_PX:
                continue

            pid += 1
            rows.append({
                "image_name": clean_image_name_from_mask(fname),
                "particle_id": pid,
                "area_pixels": area_px,
                "area_mm2": area_px * PIXEL_AREA_MM2
            })

    df = pd.DataFrame(rows)
    df.to_excel(EXCEL_AREAS, index=False)

    print("\nParticle-area table saved:", EXCEL_AREAS)


# ============================================================
# Main
# ============================================================

def main():
    run_tracer_count_and_analysis()
    run_particle_areas_from_masks()
    print("\nAll done (counts + areas).")


if __name__ == "__main__":
    main()
