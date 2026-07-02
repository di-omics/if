#!/usr/bin/env python3
"""OPTIONAL: interactive 3D review in napari.

Opens the multi-channel image, label volume, and per-lineage point
annotations. Requires a display (not part of the headless pipeline).

Usage:
    python src/review_napari.py --image data/synthetic_blastocyst.ome.tif \
                                --labels outputs/labels.tif \
                                --classified outputs/classified.tsv
"""
import argparse
import numpy as np
import pandas as pd
import tifffile


def main():
    parser = argparse.ArgumentParser(description="napari 3D review")
    parser.add_argument("--image", required=True, help="multi-channel image")
    parser.add_argument("--labels", required=True, help="label volume")
    parser.add_argument("--classified", default="outputs/classified.tsv",
                        help="classified.tsv with lineage + centroids")
    parser.add_argument("--channel-names", nargs="+",
                        default=["DAPI", "marker_A", "marker_B", "marker_C"])
    args = parser.parse_args()

    import napari

    img = tifffile.imread(str(args.image))
    labels = tifffile.imread(str(args.labels))
    df = pd.read_csv(args.classified, sep="\t")

    viewer = napari.Viewer(title="blastocyst-if review")

    # add channels
    if img.ndim == 4 and img.shape[0] <= 10:
        for c, name in enumerate(args.channel_names[:img.shape[0]]):
            viewer.add_image(img[c], name=name, blending="additive",
                             visible=(name == "DAPI"))
    else:
        viewer.add_image(img, name="image")

    # add labels
    viewer.add_labels(labels, name="nuclei")

    # add per-lineage points
    lineage_colors = {"EPI": "cyan", "PE": "green", "TE": "magenta"}
    for lin in df.lineage.unique():
        sub = df[df.lineage == lin]
        coords = sub[["centroid_z", "centroid_y", "centroid_x"]].values
        color = lineage_colors.get(lin, "white")
        viewer.add_points(coords, name=lin, face_color=color, size=6)

    napari.run()


if __name__ == "__main__":
    main()
