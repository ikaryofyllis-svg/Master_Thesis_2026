# Refactored MSc Thesis PBMC GRN Pipeline

This folder contains the cleaned/refactored version of the MSc thesis pipeline.

The original project remains outside this folder:

../

Main concepts:
- evidence-based baseline
- pseudo-supervised tabular model
- GNN-refined model
- edge-level confidence refinement
- TF-level biological prioritization
- broad motif TF artifacts

Important interpretation:
The GNN does not validate external biological ground truth.
It learns agreement with pseudo-labels derived from the evidence-based baseline.
