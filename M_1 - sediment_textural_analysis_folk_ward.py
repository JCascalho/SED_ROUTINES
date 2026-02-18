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


import numpy as np
import pandas as pd
import sys
import argparse
from typing import Dict, Tuple, List

# ---------------------------
# Unit helpers
# ---------------------------

def mm_to_phi(mm: np.ndarray) -> np.ndarray:
    return -np.log2(mm)

def phi_to_mm(phi: np.ndarray) -> np.ndarray:
    return 2 ** (-phi)

def looks_like_mm(vals: np.ndarray) -> bool:
    """Heuristic: all positive, typical range ~[0.0001, 64] mm."""
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return False
    return np.all(vals > 0) and (vals.max() <= 64) and (vals.min() >= 1e-4)

# ---------------------------
# Analyzer
# ---------------------------

class SedimentTexturalAnalyzer:
    def __init__(self):
        # Classification thresholds
        self.SORTING_THRESHOLDS = [
            (0.35, "Very well sorted"),
            (0.50, "Well sorted"),
            (0.71, "Moderately well sorted"),
            (1.00, "Moderately sorted"),
            (2.00, "Poorly sorted"),
            (4.00, "Very poorly sorted"),
            (float('inf'), "Extremely poorly sorted")
        ]
        self.SKEWNESS_THRESHOLDS = [
            (-0.30, "Very coarse skewed"),
            (-0.10, "Coarse skewed"),
            (0.10, "Symmetrical"),
            (0.30, "Fine skewed"),
            (float('inf'), "Very fine skewed")
        ]
        self.KURTOSIS_THRESHOLDS = [
            (0.67, "Very platykurtic"),
            (0.90, "Platykurtic"),
            (1.11, "Mesokurtic"),
            (1.50, "Leptokurtic"),
            (3.00, "Very leptokurtic"),
            (float('inf'), "Extremely leptokurtic")
        ]

    # ----- Percentiles & stats -----

    def calculate_percentiles(self, phi_values: np.ndarray,
                              cumulative_percentages: np.ndarray) -> Dict[str, float]:
        # Clip tail percentiles slightly to avoid end artifacts
        max_percentage = 99.98
        clipped = np.minimum(cumulative_percentages, max_percentage)
        keys = [5, 16, 25, 50, 75, 84, 95]
        return {f'phi{k}': float(np.interp(k, clipped, phi_values)) for k in keys}

    def calculate_mean_grain_size(self, phi16: float, phi50: float, phi84: float) -> float:
        return (phi16 + phi50 + phi84) / 3.0

    def calculate_sorting(self, phi5: float, phi16: float, phi84: float, phi95: float) -> float:
        return ((phi84 - phi16) / 4.0) + ((phi95 - phi5) / 6.6)

    def calculate_skewness(self, phi5: float, phi16: float, phi50: float, phi84: float, phi95: float) -> float:
        term1 = (phi84 + phi16 - 2 * phi50) / (2 * (phi84 - phi16)) if (phi84 - phi16) != 0 else np.nan
        term2 = (phi95 + phi5  - 2 * phi50) / (2 * (phi95 - phi5 )) if (phi95 - phi5 ) != 0 else np.nan
        return term1 + term2

    def calculate_kurtosis(self, phi5: float, phi25: float, phi75: float, phi95: float) -> float:
        iqr = (phi75 - phi25)
        denom = 2.44 * iqr if iqr != 0 else np.nan
        return (phi95 - phi5) / denom

    # ----- Classifications -----

    def classify_grain_size(self, mean_phi: float) -> str:
        # Bands for mean φ (Mz): >8 Clay; 4–8 Silt; 3–4 VFS; 2–3 FS; 1–2 MS; 0–1 CS; −1–0 VCS; <−1 Gravel
        if mean_phi > 8:
            return "Clay"
        elif mean_phi > 4:
            return "Silt"
        elif mean_phi > 3:
            return "Very fine sand"
        elif mean_phi > 2:
            return "Fine sand"
        elif mean_phi > 1:
            return "Medium sand"
        elif mean_phi > 0:
            return "Coarse sand"
        elif mean_phi > -1:
            return "Very coarse sand"
        else:
            return "Gravel"

    def classify_parameter(self, value: float, thresholds: List[Tuple], name: str) -> str:
        for ub, label in thresholds:
            if value <= ub:
                return label
        raise ValueError(f"Invalid {name} value: {value}")

    def classify_sorting(self, v: float) -> str:
        return self.classify_parameter(v, self.SORTING_THRESHOLDS, "sorting")

    def classify_skewness(self, v: float) -> str:
        return self.classify_parameter(v, self.SKEWNESS_THRESHOLDS, "skewness")

    def classify_kurtosis(self, v: float) -> str:
        return self.classify_parameter(v, self.KURTOSIS_THRESHOLDS, "kurtosis")

    # ----- Convenience -----

    def phi_to_millimeters(self, phi_value):
        return 2 ** (-phi_value)

    def analyze_sample(self, phi_values: np.ndarray, frequencies: pd.Series) -> Dict:
        # Robust cumulative (non-decreasing)
        freq = frequencies.astype(float).fillna(0.0).clip(lower=0.0)
        total = freq.sum()
        if total <= 0:
            # Return NaNs to be QA-flagged
            pct = {k: np.nan for k in ['phi5','phi16','phi25','phi50','phi75','phi84','phi95']}
            return {
                **pct,
                'mean_size_phi': np.nan,
                'mean_size_mm': np.nan,
                'sorting': np.nan,
                'skewness': np.nan,
                'kurtosis': np.nan,
                'grain_size_class': "N/A",
                'sorting_class': "N/A",
                'skewness_class': "N/A",
                'kurtosis_class': "N/A",
                'total_raw_percent': float(total),
                '_cum_raw': None,
                '_cum_clean': None,
            }

        cum_pct_raw = (freq.cumsum() / total) * 100.0
        cum_pct = np.maximum.accumulate(cum_pct_raw.to_numpy())
        cum_pct = np.clip(cum_pct, 0.0, 100.0)

        pct = self.calculate_percentiles(phi_values, cum_pct)

        # Stats
        phi5, phi16, phi25 = pct['phi5'], pct['phi16'], pct['phi25']
        phi50, phi75, phi84, phi95 = pct['phi50'], pct['phi75'], pct['phi84'], pct['phi95']

        if np.isnan([phi5, phi16, phi25, phi50, phi75, phi84, phi95]).any():
            mean_size_phi = sorting = skewness = kurtosis = np.nan
            classes = {
                "grain_size_class": "N/A",
                "sorting_class": "N/A",
                "skewness_class": "N/A",
                "kurtosis_class": "N/A"
            }
        else:
            mean_size_phi = self.calculate_mean_grain_size(phi16, phi50, phi84)
            sorting = self.calculate_sorting(phi5, phi16, phi84, phi95)
            skewness = self.calculate_skewness(phi5, phi16, phi50, phi84, phi95)
            kurtosis = self.calculate_kurtosis(phi5, phi25, phi75, phi95)
            classes = {
                "grain_size_class": self.classify_grain_size(mean_size_phi),
                "sorting_class":    self.classify_sorting(sorting),
                "skewness_class":   self.classify_skewness(skewness),
                "kurtosis_class":   self.classify_kurtosis(kurtosis),
            }

        return {
            **pct,
            'mean_size_phi': mean_size_phi,
            'mean_size_mm': self.phi_to_millimeters(mean_size_phi) if np.isfinite(mean_size_phi) else np.nan,
            'sorting': sorting,
            'skewness': skewness,
            'kurtosis': kurtosis,
            **classes,
            'total_raw_percent': float(total),
            '_cum_raw': cum_pct_raw.to_numpy(),
            '_cum_clean': cum_pct,
        }

# ---------------------------
# References
# ---------------------------

def create_reference_tables() -> Dict[str, pd.DataFrame]:
    phi_scale = pd.DataFrame({
        'Phi (φ)': [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8],
        'Grain Size (mm)': [32, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.031, 0.016, 0.008, 0.004],
        'Sediment Class': [
            'Cobbles', 'V coarse pebbles', 'Coarse pebbles', 'Medium pebbles',
            'Fine pebbles (granules)', 'Very coarse sand', 'Coarse sand', 'Medium sand',
            'Fine sand', 'Very fine sand', 'V coarse silt', 'Coarse silt', 'Medium silt', 'Fine silt'
        ]
    })
    classification_criteria = pd.DataFrame({
        'Parameter': ['Mean Size (φ)', 'Sorting (σI)', 'Skewness (SkI)', 'Kurtosis (KG)'],
        'Classification Ranges': [
            'Clay (>8), Silt (4–8), VFS (3–4), FS (2–3), MS (1–2), CS (0–1), VCS (−1–0), Gravel (<−1)',
            'V well (<0.35), Well (0.35–0.50), Mod well (0.50–0.71), Moderate (0.71–1.00), Poor (1.00–2.00), V poor (2.00–4.00), Ext poor (>4.00)',
            'V coarse (<−0.3), Coarse (−0.3–−0.1), Symm (−0.1–0.1), Fine (0.1–0.3), V fine (>0.3)',
            'V platy (<0.67), Platy (0.67–0.9), Meso (0.9–1.11), Lepto (1.11–1.5), V lepto (1.5–3.0), Ext lepto (>3.0)'
        ],
        'Formula': [
            '(φ16 + φ50 + φ84)/3',
            '(φ84 − φ16)/4 + (φ95 − φ5)/6.6',
            '(φ84 + φ16 − 2φ50)/(2(φ84 − φ16)) + (φ95 + φ5 − 2φ50)/(2(φ95 − φ5))',
            '(φ95 − φ5)/(2.44(φ75 − φ25))'
        ]
    })
    return {'krumbein_scale': phi_scale, 'classification_criteria': classification_criteria}

# ---------------------------
# Main
# ---------------------------

def main():
    print("=== SEDIMENT TEXTURAL ANALYSIS - FOLK & WARD METHOD (SAMPLES = ROWS) ===")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        default="INPUT_TEXTURAL_DATA_NATIVE_SEDIMENT_AND_TRACER.xlsx",
        help="Input Excel file with samples as rows and phi/mm classes as columns."
    )
    parser.add_argument(
        '-o', '--output',
        default="OUTPUT_FOLK_WARD_GRAPHICAL_METHOD_RESULTS.xlsx",
        help="Output Excel file."
    )
    parser.add_argument(
        '--units',
        choices=['auto','phi','mm'],
        default='auto',
        help="Units of grain-size classes (column headers from column 2). Default: auto-detect."
    )
    args = parser.parse_args()

    try:
        analyzer = SedimentTexturalAnalyzer()
        print(f"Loading: {args.input}")
        df = pd.read_excel(args.input)

        # ---------------------------
        # Interpret layout:
        # - first column = sample name
        # - remaining columns = grain-size classes (phi or mm) as headers
        # ---------------------------
        if df.shape[1] < 2:
            raise ValueError("Input must have at least 2 columns: sample name + one size class.")

        sample_names = df.iloc[:, 0].astype(str).tolist()
        header_classes = df.columns[1:]

        # Convert header classes to numeric sizes
        size_raw = pd.to_numeric(header_classes, errors='coerce').to_numpy()
        mask_valid = np.isfinite(size_raw)
        if not mask_valid.any():
            raise ValueError("No valid numeric grain-size classes found in column headers.")

        size_raw = size_raw[mask_valid]
        valid_cols = [c for c, ok in zip(header_classes, mask_valid) if ok]

        # Sub-dataframe with only valid class columns
        data = df[valid_cols]

        # Determine units and convert to phi
        if args.units == 'auto':
            is_mm = looks_like_mm(size_raw)
        else:
            is_mm = (args.units == 'mm')

        if is_mm:
            size_mm = size_raw
            size_phi = mm_to_phi(size_mm)
            unit_used = 'mm (converted to φ)'
        else:
            size_phi = size_raw
            unit_used = 'φ'

        # Sort by increasing phi (coarse -> fine)
        order = np.argsort(size_phi)
        phi_vals = size_phi[order]
        ordered_cols = [valid_cols[i] for i in order]
        data = data[ordered_cols]   # reorder columns

        print(f"Detected header units: {unit_used}")
        print(f"Found {len(phi_vals)} grain-size classes; {data.shape[0]} samples")

        results = []
        qa_rows = []
        percentile_cols = ['phi5','phi16','phi25','phi50','phi75','phi84','phi95']

        # Loop over rows (each row = sample)
        for sample_name, row in zip(sample_names, data.itertuples(index=False)):
            freqs = pd.to_numeric(pd.Series(row, index=ordered_cols), errors='coerce').fillna(0.0)
            total = float(freqs.sum())

            # Analyze (returns percentiles & stats + raw/clean cumulative for QA)
            out = analyzer.analyze_sample(phi_vals, freqs)

            # Build results row
            results.append({
                'sample': sample_name,
                **{k: out.get(k, np.nan) for k in percentile_cols},
                'mean_size_phi': out['mean_size_phi'],
                'mean_size_mm': out['mean_size_mm'],
                'sorting': out['sorting'],
                'skewness': out['skewness'],
                'kurtosis': out['kurtosis'],
                'grain_size_class': out['grain_size_class'],
                'sorting_class': out['sorting_class'],
                'skewness_class': out['skewness_class'],
                'kurtosis_class': out['kurtosis_class'],
            })

            # ---------------------------
            # QA checks
            # ---------------------------
            issues = []
            has_nan_pcts = np.isnan([out.get(c) for c in percentile_cols]).any()

            # Raw cumulative monotonicity
            if total > 0 and isinstance(out.get('_cum_raw'), np.ndarray):
                cum_raw = out['_cum_raw']
                diffs = np.diff(cum_raw[~np.isnan(cum_raw)])
                raw_monotonic = bool(np.all(diffs >= -1e-8))
                if not raw_monotonic:
                    issues.append("Raw cumulative not monotonic (check input ordering).")
            else:
                raw_monotonic = False
                issues.append("No material (sum = 0) or missing cumulative curve.")

            # Percentile-based checks
            phi5  = out.get('phi5')
            phi16 = out.get('phi16')
            phi25 = out.get('phi25')
            phi50 = out.get('phi50')   # φ50 now used explicitly in QA
            phi75 = out.get('phi75')
            phi84 = out.get('phi84')
            phi95 = out.get('phi95')

            center_zero = (np.isfinite(phi84) and np.isfinite(phi16) and (phi84 - phi16) == 0)
            tail_zero   = (np.isfinite(phi95) and np.isfinite(phi5)  and (phi95 - phi5)   == 0)
            iqr_zero    = (np.isfinite(phi75) and np.isfinite(phi25) and (phi75 - phi25) == 0)

            if center_zero:
                issues.append("φ84 − φ16 = 0 (no central spread).")
            if tail_zero:
                issues.append("φ95 − φ5 = 0 (no tail spread).")
            if iqr_zero:
                issues.append("φ75 − φ25 = 0 (IQR = 0; KG undefined).")
            if has_nan_pcts:
                issues.append("One or more percentiles are NaN (insufficient spread or empty bins).")

            # φ50 consistency
            if not np.isfinite(phi50):
                issues.append("φ50 (median) undefined or non-finite.")
            else:
                if np.isfinite(phi16) and phi50 < phi16:
                    issues.append("φ50 (median) < φ16 (check percentile curve).")
                if np.isfinite(phi84) and phi50 > phi84:
                    issues.append("φ50 (median) > φ84 (check percentile curve).")

            qa_rows.append({
                'Sample': sample_name,
                'Units_Used': unit_used,
                'Total_raw_sum_%': total,
                'Raw_Cumulative_Monotonic': raw_monotonic,
                'Has_NaN_percentiles': bool(has_nan_pcts),
                'Center_width_zero': bool(center_zero),
                'Tail_width_zero': bool(tail_zero),
                'IQR_zero': bool(iqr_zero),
                'Issues': "; ".join(issues) if issues else "OK"
            })

        # Assemble results DataFrame
        results_df = pd.DataFrame(results)

        # Also export percentiles in mm (for convenience)
        mm_df = results_df[['sample']].copy()
        for c in percentile_cols:
            if c in results_df.columns:
                mm_df[c.replace('phi','d_mm_at_')] = phi_to_mm(results_df[c])

        # Summary rows (means)
        mean_phi_row = {
            'sample': 'MEAN',
            **{c: results_df[c].mean() for c in results_df.select_dtypes(np.number).columns}
        }
        mean_phi_df = pd.DataFrame([mean_phi_row])

        mean_mm_row = {
            'sample': 'MEAN',
            **{c: mm_df[c].mean() for c in mm_df.columns if c != 'sample'}
        }
        mean_mm_df = pd.DataFrame([mean_mm_row])

        qa_df = pd.DataFrame(qa_rows)

        # References
        refs = create_reference_tables()

        # Write Excel
        with pd.ExcelWriter(args.output, engine='openpyxl') as w:
            results_df.to_excel(w, sheet_name='Textural_Parameters_phi', index=False)
            mean_phi_df.to_excel(
                w,
                sheet_name='Textural_Parameters_phi',
                startrow=len(results_df)+1,
                header=False,
                index=False
            )

            mm_df.to_excel(w, sheet_name='Percentiles_mm', index=False)
            mean_mm_df.to_excel(
                w,
                sheet_name='Percentiles_mm',
                startrow=len(mm_df)+1,
                header=False,
                index=False
            )

            qa_df.to_excel(w, sheet_name='QA', index=False)

            refs['krumbein_scale'].to_excel(w, sheet_name='Phi_Scale_Reference', index=False)
            refs['classification_criteria'].to_excel(w, sheet_name='Classification_Criteria', index=False)

        print("OK — saved:", args.output)

    except Exception as e:
        print("ERROR:", str(e))
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
