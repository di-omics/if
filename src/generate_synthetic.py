#!/usr/bin/env python3
"""Generate a synthetic 3D multi-channel OME-TIFF for pipeline testing.

Creates a DAPI-like nuclear channel plus 3 marker channels (generic names:
marker_A, marker_B, marker_C) with embedded nuclei of varying intensity.
Also writes a ground-truth label volume so downstream steps can be tested
without Cellpose / GPU.

All data is synthetic.
"""
import numpy as np
import tifffile
from pathlib import Path
from scipy.ndimage import gaussian_filter

RNG = np.random.default_rng(42)
OUT = Path(__file__).resolve().parents[1] / "data"
OUT.mkdir(exist_ok=True)

# Volume dimensions (Z, Y, X)
NZ, NY, NX = 30, 256, 256
N_NUCLEI = 40
RADIUS_RANGE = (8, 16)  # nucleus radius in pixels (XY)
Z_RADIUS_SCALE = 0.4    # nuclei are flatter in Z

CHANNELS = ["DAPI", "marker_A", "marker_B", "marker_C"]


def make_sphere_mask(shape, center, radius_yx, radius_z):
    """Binary mask for an ellipsoid nucleus."""
    zz, yy, xx = np.ogrid[:shape[0], :shape[1], :shape[2]]
    cz, cy, cx = center
    dist = ((zz - cz) / radius_z) ** 2 + ((yy - cy) / radius_yx) ** 2 + ((xx - cx) / radius_yx) ** 2
    return dist <= 1.0


def main():
    volume = np.zeros((len(CHANNELS), NZ, NY, NX), dtype=np.float32)
    labels = np.zeros((NZ, NY, NX), dtype=np.int32)

    # assign each nucleus a "lineage" for marker intensities
    # roughly: ~40% TE-like, ~30% EPI-like, ~30% PE-like
    lineage_probs = [0.4, 0.3, 0.3]
    lineage_names = ["TE", "EPI", "PE"]

    nuclei_info = []
    for i in range(1, N_NUCLEI + 1):
        r_yx = RNG.integers(RADIUS_RANGE[0], RADIUS_RANGE[1])
        r_z = max(2, int(r_yx * Z_RADIUS_SCALE))
        cy = RNG.integers(r_yx + 2, NY - r_yx - 2)
        cx = RNG.integers(r_yx + 2, NX - r_yx - 2)
        cz = RNG.integers(r_z + 1, NZ - r_z - 1)

        mask = make_sphere_mask((NZ, NY, NX), (cz, cy, cx), r_yx, r_z)
        # avoid overlapping labels (first writer wins)
        overlap = labels[mask] > 0
        if overlap.any():
            mask_flat = mask.copy()
            mask_flat[labels > 0] = False
            if mask_flat.sum() < 20:
                continue
            mask = mask_flat

        labels[mask] = i

        # DAPI: all nuclei bright
        dapi_intensity = RNG.uniform(0.6, 1.0)
        volume[0][mask] += dapi_intensity

        # marker intensities depend on lineage
        lineage = RNG.choice(len(lineage_names), p=lineage_probs)
        lname = lineage_names[lineage]

        # marker_A high in EPI, low elsewhere
        # marker_B high in PE, low elsewhere
        # marker_C high in TE, low elsewhere
        marker_levels = {
            "EPI": [0.8, 0.15, 0.1],
            "PE":  [0.1, 0.8, 0.15],
            "TE":  [0.1, 0.1, 0.75],
        }
        for ch_idx, level in enumerate(marker_levels[lname]):
            intensity = level * RNG.uniform(0.7, 1.3)
            volume[ch_idx + 1][mask] += intensity

        nuclei_info.append({
            "label": i, "lineage": lname,
            "z": cz, "y": cy, "x": cx, "r_yx": r_yx,
        })

    # add Gaussian blur + Poisson-like noise for realism
    for c in range(len(CHANNELS)):
        volume[c] = gaussian_filter(volume[c], sigma=(1.0, 1.5, 1.5))
        noise = RNG.normal(0, 0.02, volume[c].shape).astype(np.float32)
        volume[c] = np.clip(volume[c] + noise, 0, None)

    # scale to uint16 range
    for c in range(len(CHANNELS)):
        mx = volume[c].max()
        if mx > 0:
            volume[c] = volume[c] / mx
    img_u16 = (volume * 65535).astype(np.uint16)

    # write OME-TIFF (CZYX)
    out_img = OUT / "synthetic_blastocyst.ome.tif"
    ome_metadata = {
        "axes": "CZYX",
        "Channel": {"Name": CHANNELS},
    }
    tifffile.imwrite(
        str(out_img), img_u16, ome=True, photometric="minisblack",
        metadata=ome_metadata,
    )

    # write ground-truth labels
    out_labels = OUT / "synthetic_labels.tif"
    tifffile.imwrite(str(out_labels), labels)

    # write ground-truth lineage table
    import pandas as pd
    pd.DataFrame(nuclei_info).to_csv(OUT / "synthetic_truth.tsv", sep="\t", index=False)

    n_actual = len(nuclei_info)
    print(f"wrote {out_img}  ({len(CHANNELS)}C x {NZ}Z x {NY}Y x {NX}X, {n_actual} nuclei)")
    print(f"wrote {out_labels}  (label volume, {n_actual} labels)")
    print(f"wrote {OUT / 'synthetic_truth.tsv'}  (ground-truth lineage)")


if __name__ == "__main__":
    main()
