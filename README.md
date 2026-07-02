# blastocyst-if

3D immunofluorescence analysis pipeline for whole-mount preimplantation
blastocysts: Cellpose nuclear segmentation, per-nucleus marker measurement,
and generic lineage classification.

All default data is synthetic. The pipeline runs end-to-end with no GPU
and no data download.

## Stack

- **Cellpose / Cellpose-SAM** - 3D nuclear segmentation (MPS on Apple Silicon)
- **scikit-image** - per-nucleus regionprops (volume, centroid, intensity)
- **napari** - optional interactive 3D review
- **bioio** - CZI / LIF / OME-TIFF readers

## Setup

```bash
conda env create -f environment.yml
conda activate blastocyst-if
```

## One-command synthetic run

```bash
make synthetic
```

This runs the full chain on a synthetic 3D multi-channel OME-TIFF:

1. **generate_synthetic.py** - writes a DAPI + 3 marker channel volume with
   embedded nuclei of varying intensity, plus a ground-truth label volume.
2. **01_segment.py --synthetic** - loads the ground-truth labels (skips
   Cellpose so no GPU is needed).
3. **02_measure.py** - per-nucleus volume, centroid, and per-channel
   mean/integrated intensity via `skimage.regionprops`.
4. **03_classify.py** - assigns each nucleus a lineage (EPI / PE / TE)
   based on which marker channel has the highest normalized intensity.
   A `--method gmm` flag switches to a 2-component GMM gate.
5. **plots.py** - static QC figure: nucleus volume distribution, marker
   intensities by lineage, and classification proportions.

### Optional: napari review

```bash
python src/review_napari.py \
    --image data/synthetic_blastocyst.ome.tif \
    --labels outputs/labels.tif \
    --classified outputs/classified.tsv
```

Opens the image, labels, and per-lineage centroid points in napari.
Requires a display.

## Run on real data

Point `01_segment.py` at a downloaded CZI, LIF, or OME-TIFF and specify
the DAPI channel index:

```bash
python src/01_segment.py --input /path/to/blastocyst.czi --dapi-channel 0
python src/02_measure.py --image /path/to/blastocyst.czi --labels outputs/labels.tif \
    --channel-names DAPI NANOG GATA6 CDX2
python src/03_classify.py --input outputs/measurements.tsv \
    --markers NANOG_mean GATA6_mean CDX2_mean \
    --lineages EPI PE TE
python src/plots.py
```

Published confocal stacks are available from:
- **Niakan lab, Simon et al. 2025** (Nat Commun) - figshare [28597145](https://doi.org/10.6084/m9.figshare.28597145)
- Reference pipeline - zenodo [15640446](https://doi.org/10.5281/zenodo.15640446)

![QC plots](assets/blastocyst_qc.png)

## Files

```
src/generate_synthetic.py   synthetic 3D multi-channel OME-TIFF + ground-truth labels
src/01_segment.py           Cellpose 3D nuclear segmentation (or --synthetic for GT)
src/02_measure.py           per-nucleus measurements via regionprops
src/03_classify.py          generic marker-threshold / GMM lineage classification
src/plots.py                static QC figure
src/review_napari.py        optional interactive 3D review in napari
```
