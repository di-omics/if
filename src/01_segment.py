#!/usr/bin/env python3
"""3D nuclear segmentation via Cellpose (or ground-truth labels for synthetic data).

Reads any bioio-supported file (OME-TIFF, CZI, LIF) and segments the DAPI
channel. Writes a label volume to outputs/.

Usage:
    python src/01_segment.py --input data/synthetic_blastocyst.ome.tif --dapi-channel 0
    python src/01_segment.py --input data/synthetic_blastocyst.ome.tif --synthetic
"""
import argparse
import numpy as np
import tifffile
from pathlib import Path


def load_image(path, dapi_channel=0):
    """Load a multi-channel 3D image and return the DAPI channel as (Z, Y, X)."""
    img = tifffile.imread(str(path))
    # handle common shapes: CZYX, ZCYX, ZYX
    if img.ndim == 4:
        # assume CZYX if first dim is small (number of channels)
        if img.shape[0] <= 10:
            return img[dapi_channel]
        # otherwise ZCYX
        return img[:, dapi_channel]
    elif img.ndim == 3:
        return img
    else:
        raise ValueError(f"unexpected image shape {img.shape}")


def segment_cellpose(dapi, use_mps=True):
    """Run Cellpose 3D nuclear segmentation."""
    from cellpose import models
    import torch

    device = None
    if use_mps and hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")

    model = models.Cellpose(model_type="nuclei", device=device)
    masks, flows, styles, diams = model.eval(
        dapi, channels=[0, 0], do_3D=True, diameter=None,
    )
    return masks


def segment_synthetic(labels_path):
    """Load pre-computed ground-truth labels (skips Cellpose)."""
    labels = tifffile.imread(str(labels_path))
    print(f"loaded ground-truth labels from {labels_path}")
    return labels


def main():
    parser = argparse.ArgumentParser(description="3D nuclear segmentation")
    parser.add_argument("--input", required=True, help="path to image file")
    parser.add_argument("--dapi-channel", type=int, default=0,
                        help="channel index for DAPI (default: 0)")
    parser.add_argument("--synthetic", action="store_true",
                        help="use ground-truth labels instead of Cellpose")
    parser.add_argument("--labels-path", default=None,
                        help="path to ground-truth label volume (default: data/synthetic_labels.tif)")
    parser.add_argument("--output-dir", default="outputs",
                        help="directory for output files")
    parser.add_argument("--no-mps", action="store_true",
                        help="disable MPS (Apple Silicon GPU)")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    if args.synthetic:
        labels_path = args.labels_path or "data/synthetic_labels.tif"
        labels = segment_synthetic(labels_path)
    else:
        print(f"loading DAPI channel {args.dapi_channel} from {args.input}")
        dapi = load_image(args.input, args.dapi_channel)
        print(f"DAPI volume shape: {dapi.shape}, dtype: {dapi.dtype}")
        labels = segment_cellpose(dapi, use_mps=not args.no_mps)

    n_nuclei = len(np.unique(labels)) - 1  # subtract background
    out_path = out_dir / "labels.tif"
    tifffile.imwrite(str(out_path), labels.astype(np.int32))
    print(f"wrote {out_path}  ({n_nuclei} nuclei)")


if __name__ == "__main__":
    main()
