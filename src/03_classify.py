#!/usr/bin/env python3
"""Generic marker-threshold lineage classification.

Assigns each nucleus a lineage label (e.g. EPI / PE / TE) based on which
marker channel has the highest normalized mean intensity. Optionally fits a
2-component GMM per marker to set an automatic high/low gate.

This is a GENERIC method -- no hardcoded thresholds, channel maps, or
embryo-specific logic. It works on any measurements.tsv that has per-channel
mean intensity columns.

Usage:
    python src/03_classify.py --input outputs/measurements.tsv
    python src/03_classify.py --input outputs/measurements.tsv --method gmm
"""
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def classify_max_marker(df, marker_cols, lineage_names):
    """Assign lineage by highest normalized marker intensity."""
    marker_vals = df[marker_cols].values.astype(float)
    # normalize each marker to [0, 1] across all nuclei
    mins = marker_vals.min(axis=0)
    maxs = marker_vals.max(axis=0)
    ranges = maxs - mins
    ranges[ranges == 0] = 1
    normed = (marker_vals - mins) / ranges
    best = np.argmax(normed, axis=1)
    return [lineage_names[i] for i in best]


def classify_gmm(df, marker_cols, lineage_names):
    """Assign lineage via 2-component GMM per marker, then take highest posterior."""
    from sklearn.mixture import GaussianMixture

    posteriors = np.zeros((len(df), len(marker_cols)))
    for i, col in enumerate(marker_cols):
        vals = df[col].values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=0).fit(vals)
        # the component with the higher mean is the "positive" component
        high_comp = np.argmax(gmm.means_.ravel())
        posteriors[:, i] = gmm.predict_proba(vals)[:, high_comp]
    best = np.argmax(posteriors, axis=1)
    return [lineage_names[i] for i in best]


def main():
    parser = argparse.ArgumentParser(description="lineage classification")
    parser.add_argument("--input", required=True, help="measurements.tsv")
    parser.add_argument("--markers", nargs="+",
                        default=["marker_A_mean", "marker_B_mean", "marker_C_mean"],
                        help="marker columns to use for classification")
    parser.add_argument("--lineages", nargs="+", default=["lineage_A", "lineage_B", "lineage_C"],
                        help="lineage names (same order as markers)")
    parser.add_argument("--method", choices=["max", "gmm"], default="max",
                        help="classification method (default: max)")
    parser.add_argument("--output-dir", default="outputs")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    marker_cols = args.markers
    lineage_names = args.lineages

    if args.method == "gmm":
        df["lineage"] = classify_gmm(df, marker_cols, lineage_names)
    else:
        df["lineage"] = classify_max_marker(df, marker_cols, lineage_names)

    # proportions
    counts = df["lineage"].value_counts()
    total = len(df)

    print(f"\n=== Lineage classification ({args.method} method) ===")
    print(f"nuclei: {total}")
    for lin in lineage_names:
        n = counts.get(lin, 0)
        print(f"  {lin}: {n} ({100 * n / total:.1f}%)")

    out_classified = out_dir / "classified.tsv"
    df.to_csv(out_classified, sep="\t", index=False)

    props = pd.DataFrame({
        "lineage": lineage_names,
        "count": [counts.get(ln, 0) for ln in lineage_names],
        "fraction": [counts.get(ln, 0) / total for ln in lineage_names],
    })
    out_props = out_dir / "proportions.tsv"
    props.to_csv(out_props, sep="\t", index=False)

    print(f"\nwrote {out_classified}, {out_props}")


if __name__ == "__main__":
    main()
