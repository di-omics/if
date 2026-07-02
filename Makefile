PY ?= python

.PHONY: synthetic clean

synthetic:
	$(PY) src/generate_synthetic.py
	$(PY) src/01_segment.py --input data/synthetic_blastocyst.ome.tif --synthetic
	$(PY) src/02_measure.py --image data/synthetic_blastocyst.ome.tif --labels outputs/labels.tif
	$(PY) src/03_classify.py --input outputs/measurements.tsv

clean:
	rm -rf data/ outputs/ assets/ __pycache__/ src/__pycache__/
