
# -*- coding: utf-8 -*-
"""
Step 04C — scATAC gene activity construction.

Purpose:
- Load clustered ATAC object from Step 04B.
- Load peak-to-gene annotation table.
- Match ATAC peak IDs to annotation peak IDs.
- Build a sparse peak × gene mapping matrix.
- Compute cells × genes gene activity matrix.
- Normalize/log1p gene activity for marker-based interpretation.
- Save gene activity AnnData, audit tables, marker plots, and manifest.

Important:
- This step does not replace RNA annotation.
- Gene activity is used to interpret ATAC clusters and support RNA–ATAC matching.
- All outputs must stay inside the active project.
"""

from pathlib import Path
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from scipy import sparse
from anndata import AnnData


try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

src_dir = refactor_dir / "src"  # reusable package folder

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # allow local pbmcgrn imports

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


def assert_path_inside_project(path: Path, project_root: Path, label: str) -> None:
    """Fail if a path is outside active project."""

    resolved_path = Path(path).resolve()  # normalize target path
    resolved_project_root = Path(project_root).resolve()  # normalize project root

    try:
        resolved_path.relative_to(resolved_project_root)  # validate containment
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside active project.\n"
            f"Path: {resolved_path}\n"
            f"Project root: {resolved_project_root}"
        ) from exc


def resolve_project_path(paths: ProjectPaths, relative_path: str, label: str) -> Path:
    """Resolve a YAML relative path inside the active project."""

    out_path = paths.prepare_output_path(relative_path, label=label)  # project-aware path

    assert_path_inside_project(out_path, paths.project_root, label)  # safety check

    return out_path


def save_show_close(path: Path, dpi: int = 300, show: bool = True) -> None:
    """Save current matplotlib figure and optionally show it."""

    plt.savefig(path, dpi=dpi, bbox_inches="tight")  # save current figure

    if show:
        plt.show()  # display interactively

    plt.close()  # close figure to avoid memory buildup


def detect_column(columns, candidates, label: str):
    """Detect a column from a list of candidate names."""

    lower_to_original = {col.lower(): col for col in columns}  # case-insensitive lookup

    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]  # return original column name

    raise ValueError(
        f"Could not detect {label} column.\n"
        f"Available columns: {list(columns)}\n"
        f"Tried candidates: {candidates}"
    )


def normalize_peak_id(value: str) -> str:
    """Normalize peak IDs to a common chr:start-end style where possible."""

    value = str(value).strip()  # remove whitespace

    value = value.replace("_", ":") if value.count("_") == 1 and ":" not in value else value  # light normalization

    return value


def build_peak_gene_mapping(atac_peak_ids, annotation, peak_col, gene_col, min_peaks_per_gene: int):
    """Build sparse peak × gene mapping matrix from peak annotation."""

    peak_index = pd.Series(
        np.arange(len(atac_peak_ids)),
        index=pd.Index([normalize_peak_id(p) for p in atac_peak_ids], name="peak_id"),
    )  # map normalized ATAC peak ID to column index

    annot = annotation[[peak_col, gene_col]].copy()  # keep only required columns

    annot[peak_col] = annot[peak_col].map(normalize_peak_id)  # normalize annotation peak IDs

    annot[gene_col] = annot[gene_col].astype(str).str.strip()  # clean gene names

    annot = annot.dropna(subset=[peak_col, gene_col])  # remove missing mapping rows

    annot = annot[annot[gene_col] != ""]  # remove empty gene names

    annot = annot[annot[gene_col].str.lower() != "nan"]  # remove string nan values

    annot["peak_idx"] = annot[peak_col].map(peak_index)  # match annotation peaks to ATAC var_names

    matched = annot.dropna(subset=["peak_idx"]).copy()  # keep only matched peaks

    matched["peak_idx"] = matched["peak_idx"].astype(int)  # peak matrix column index

    matched = matched.drop_duplicates(subset=["peak_idx", gene_col])  # avoid duplicate peak-gene edges

    gene_counts = matched[gene_col].value_counts()  # count peaks per gene

    keep_genes = gene_counts[gene_counts >= min_peaks_per_gene].index  # genes with enough mapped peaks

    matched = matched[matched[gene_col].isin(keep_genes)].copy()  # filter weakly mapped genes

    genes = pd.Index(sorted(matched[gene_col].unique()), name="gene")  # stable gene order

    gene_index = pd.Series(np.arange(len(genes)), index=genes)  # gene name to matrix column index

    row_idx = matched["peak_idx"].to_numpy()  # peak indices

    col_idx = matched[gene_col].map(gene_index).to_numpy()  # gene indices

    data = np.ones(len(matched), dtype=np.float32)  # binary peak-to-gene edge weights

    mapping = sparse.csr_matrix(
        (data, (row_idx, col_idx)),
        shape=(len(atac_peak_ids), len(genes)),
    )  # sparse peak × gene matrix

    return mapping, genes, matched


# ============================================================
# 0. Project setup
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # create standard output dirs

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 1. Load config
# ============================================================

config_path = paths.configs / "gene_activity_atac.yaml"  # Step 04C config

assert_path_inside_project(config_path, paths.project_root, "gene_activity_atac.yaml")  # safety check

if not config_path.exists():
    raise FileNotFoundError(f"gene_activity_atac.yaml not found:\n{config_path}")

config = load_yaml(config_path)  # load selected config

input_cfg = config["input"]  # input paths
annotation_cfg = config["annotation"]  # annotation column settings
gene_activity_cfg = config["gene_activity"]  # gene activity settings
markers_cfg = config["markers"]  # marker genes
outputs_cfg = config["outputs"]  # output paths

print("\nLoaded gene activity config:")
print(config_path)


# ============================================================
# 2. Resolve paths
# ============================================================

atac_in = paths.project_root / input_cfg["atac_clustered_h5ad"]  # clustered ATAC object

peak_annotation_in = paths.project_root / input_cfg["peak_annotation_tsv"]  # peak annotation table

gene_activity_out = resolve_project_path(
    paths,
    outputs_cfg["gene_activity_h5ad"],
    "gene_activity_h5ad",
)  # output gene activity object

mapping_summary_out = resolve_project_path(
    paths,
    outputs_cfg["mapping_summary"],
    "mapping_summary",
)  # mapping audit table

marker_availability_out = resolve_project_path(
    paths,
    outputs_cfg["marker_availability"],
    "marker_availability",
)  # marker availability table

marker_mean_by_cluster_out = resolve_project_path(
    paths,
    outputs_cfg["marker_mean_by_cluster"],
    "marker_mean_by_cluster",
)  # marker mean table

marker_dotplot_out = resolve_project_path(
    paths,
    outputs_cfg["marker_dotplot"],
    "marker_dotplot",
)  # marker dotplot figure

marker_umap_out = resolve_project_path(
    paths,
    outputs_cfg["marker_umap"],
    "marker_umap",
)  # marker UMAP figure

manifest_out = resolve_project_path(
    paths,
    outputs_cfg["manifest"],
    "step04C_manifest",
)  # manifest path

for label, path in {
    "atac_in": atac_in,
    "peak_annotation_in": peak_annotation_in,
    "gene_activity_out": gene_activity_out,
    "mapping_summary_out": mapping_summary_out,
    "marker_availability_out": marker_availability_out,
    "marker_mean_by_cluster_out": marker_mean_by_cluster_out,
    "marker_dotplot_out": marker_dotplot_out,
    "marker_umap_out": marker_umap_out,
    "manifest_out": manifest_out,
}.items():
    assert_path_inside_project(path, paths.project_root, label)  # keep project-local

gene_activity_out.parent.mkdir(parents=True, exist_ok=True)  # create output folder
mapping_summary_out.parent.mkdir(parents=True, exist_ok=True)  # create table folder
marker_dotplot_out.parent.mkdir(parents=True, exist_ok=True)  # create figure folder
manifest_out.parent.mkdir(parents=True, exist_ok=True)  # create manifest folder

if not atac_in.exists():
    raise FileNotFoundError(f"Clustered ATAC object not found:\n{atac_in}")

if not peak_annotation_in.exists():
    raise FileNotFoundError(f"Peak annotation table not found:\n{peak_annotation_in}")


# ============================================================
# 3. Load clustered ATAC and peak annotation
# ============================================================

atac = sc.read_h5ad(atac_in)  # load Step 04B clustered ATAC object

print("\nLoaded clustered ATAC object:")
print(atac)

if "leiden_0_6" not in atac.obs.columns:
    raise KeyError("leiden_0_6 not found in ATAC obs. Check Step 04B output.")

if "X_umap" not in atac.obsm.keys():
    raise KeyError("X_umap not found in ATAC obsm. Check Step 04B output.")

peak_annotation = pd.read_csv(peak_annotation_in, sep="\t")  # load peak annotation TSV

print("\nLoaded peak annotation:")
print(peak_annotation.shape)
print(peak_annotation.columns.tolist())


# ============================================================
# 4. Detect annotation columns
# ============================================================

peak_col = annotation_cfg.get("peak_column")  # optional explicit peak column from YAML

gene_col = annotation_cfg.get("gene_column")  # optional explicit gene column from YAML

build_peak_from_columns = bool(annotation_cfg.get("build_peak_from_columns", False))  # whether to construct peak IDs

if peak_col is None and build_peak_from_columns:
    chrom_col = annotation_cfg.get("chrom_column", "chrom")  # chromosome column
    start_col = annotation_cfg.get("start_column", "start")  # peak start column
    end_col = annotation_cfg.get("end_column", "end")  # peak end column

    for required_col in [chrom_col, start_col, end_col]:
        if required_col not in peak_annotation.columns:
            raise KeyError(f"Required peak coordinate column not found: {required_col}")

    peak_annotation["peak_id"] = (
        peak_annotation[chrom_col].astype(str)
        + ":"
        + peak_annotation[start_col].astype(int).astype(str)
        + "-"
        + peak_annotation[end_col].astype(int).astype(str)
    )  # build chr:start-end peak IDs to match ATAC var_names

    peak_col = "peak_id"  # use constructed peak ID column

elif peak_col is None:
    peak_col = detect_column(
        peak_annotation.columns,
        candidates=["peak", "peaks", "peak_id", "peak_name", "Peak", "name"],
        label="peak",
    )  # detect existing peak column

if gene_col is None:
    gene_col = detect_column(
        peak_annotation.columns,
        candidates=["gene", "gene_name", "symbol", "gene_symbol", "Gene", "name2"],
        label="gene",
    )  # detect gene column

if peak_col not in peak_annotation.columns:
    raise KeyError(f"Configured peak_column not found: {peak_col}")

if gene_col not in peak_annotation.columns:
    raise KeyError(f"Configured gene_column not found: {gene_col}")

print("\nUsing annotation columns:")
print("peak_col:", peak_col)
print("gene_col:", gene_col)


# ============================================================
# 5. Build peak × gene mapping
# ============================================================

min_peaks_per_gene = int(gene_activity_cfg.get("min_peaks_per_gene", 1))  # mapping filter

mapping, genes, matched_mapping = build_peak_gene_mapping(
    atac_peak_ids=atac.var_names,
    annotation=peak_annotation,
    peak_col=peak_col,
    gene_col=gene_col,
    min_peaks_per_gene=min_peaks_per_gene,
)  # construct sparse peak × gene mapping

n_atac_peaks = int(atac.n_vars)  # peaks in ATAC object
n_annotation_rows = int(peak_annotation.shape[0])  # rows in annotation table
n_matched_peak_gene_edges = int(matched_mapping.shape[0])  # matched peak-gene pairs
n_matched_peaks = int(matched_mapping["peak_idx"].nunique())  # matched unique peaks
n_genes = int(len(genes))  # output genes

if n_matched_peak_gene_edges == 0:
    raise ValueError(
        "No ATAC peaks matched the peak annotation table. "
        "Check peak ID format and peak_column."
    )

mapping_summary = pd.DataFrame(
    {
        "metric": [
            "n_atac_peaks",
            "n_annotation_rows",
            "n_matched_peak_gene_edges",
            "n_matched_unique_peaks",
            "fraction_atac_peaks_matched",
            "n_gene_activity_genes",
            "min_peaks_per_gene",
            "peak_column",
            "gene_column",
        ],
        "value": [
            n_atac_peaks,
            n_annotation_rows,
            n_matched_peak_gene_edges,
            n_matched_peaks,
            n_matched_peaks / n_atac_peaks,
            n_genes,
            min_peaks_per_gene,
            peak_col,
            gene_col,
        ],
    }
)  # mapping audit summary

write_csv(mapping_summary, mapping_summary_out, index=False)  # save mapping summary

print("\nGene activity mapping summary:")
print(mapping_summary)


# ============================================================
# 6. Compute cells × genes gene activity matrix
# ============================================================

use_layer = gene_activity_cfg.get("use_layer", "counts")  # counts layer or X

if use_layer is not None and use_layer in atac.layers.keys():
    X_peaks = atac.layers[use_layer]  # use preserved raw counts layer
else:
    X_peaks = atac.X  # fallback to current X

if not sparse.issparse(X_peaks):
    X_peaks = sparse.csr_matrix(X_peaks)  # ensure sparse multiplication
else:
    X_peaks = X_peaks.tocsr()  # efficient sparse format

X_gene_activity = X_peaks @ mapping  # cells × peaks times peaks × genes = cells × genes

X_gene_activity = X_gene_activity.tocsr()  # ensure CSR matrix for AnnData

gene_activity = AnnData(
    X=X_gene_activity,
    obs=atac.obs.copy(),
    var=pd.DataFrame(index=genes),
)  # create gene activity AnnData

gene_activity.obsm["X_umap"] = atac.obsm["X_umap"].copy()  # copy ATAC UMAP for visualization

gene_activity.uns["source"] = {
    "source_atac_h5ad": str(atac_in),
    "peak_annotation_tsv": str(peak_annotation_in),
    "peak_column": peak_col,
    "gene_column": gene_col,
    "method": "peak-to-gene summation using annotation-derived sparse mapping",
}  # provenance metadata

gene_activity.layers["raw_gene_activity_counts"] = gene_activity.X.copy()  # preserve raw activity counts

if bool(gene_activity_cfg.get("normalize_total", True)):
    sc.pp.normalize_total(
        gene_activity,
        target_sum=float(gene_activity_cfg.get("target_sum", 10000.0)),
    )  # normalize per cell for comparable marker activity

if bool(gene_activity_cfg.get("log1p", True)):
    sc.pp.log1p(gene_activity)  # log-transform for visualization and marker summaries


# ============================================================
# 7. Marker availability and marker summaries
# ============================================================

marker_genes = list(markers_cfg.get("genes", []))  # marker genes from YAML

available_markers = [gene for gene in marker_genes if gene in gene_activity.var_names]  # present markers

missing_markers = [gene for gene in marker_genes if gene not in gene_activity.var_names]  # absent markers

marker_availability = pd.DataFrame(
    {
        "marker_gene": marker_genes,
        "available_in_gene_activity": [gene in gene_activity.var_names for gene in marker_genes],
    }
)  # marker availability audit table

write_csv(marker_availability, marker_availability_out, index=False)  # save marker availability

print("\nMarker availability:")
print(marker_availability)

if len(available_markers) > 0:
    marker_df = sc.get.obs_df(
        gene_activity,
        keys=["leiden_0_6"] + available_markers,
    )  # extract marker activity and cluster labels

    marker_mean_by_cluster = (
        marker_df
        .groupby("leiden_0_6", observed=True)[available_markers]
        .mean()
        .reset_index()
    )  # mean marker activity per ATAC cluster

    write_csv(marker_mean_by_cluster, marker_mean_by_cluster_out, index=False)  # save summary table

    print("\nMarker mean gene activity by ATAC cluster:")
    print(marker_mean_by_cluster)

    sc.pl.dotplot(
        gene_activity,
        var_names=available_markers,
        groupby="leiden_0_6",
        standard_scale="var",
        show=False,
    )  # marker gene activity dotplot by ATAC clusters

    save_show_close(marker_dotplot_out, dpi=300, show=True)  # save dotplot

    umap_markers = available_markers[:12]  # avoid giant UMAP panel

    sc.pl.umap(
        gene_activity,
        color=["leiden_0_6"] + umap_markers,
        ncols=4,
        cmap="viridis",
        frameon=False,
        show=False,
    )  # UMAP of ATAC clusters and selected marker gene activities

    save_show_close(marker_umap_out, dpi=300, show=True)  # save UMAP marker figure

else:
    marker_mean_by_cluster = pd.DataFrame()  # empty fallback if no markers found

    write_csv(marker_mean_by_cluster, marker_mean_by_cluster_out, index=False)  # save empty table

    print("\nWarning: No requested marker genes were found in gene activity var_names.")


# ============================================================
# 8. Save gene activity object
# ============================================================

gene_activity.write_h5ad(gene_activity_out)  # save final gene activity object

print("\nSaved gene activity object to:")
print(gene_activity_out)


# ============================================================
# 9. Manifest
# ============================================================

manifest = {
    "step": "04C",
    "description": "scATAC gene activity construction from peak accessibility and peak-to-gene annotation",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "clustered_atac_h5ad": str(atac_in),
        "peak_annotation_tsv": str(peak_annotation_in),
        "gene_activity_config": str(config_path),
    },
    "annotation_columns": {
        "peak_column": peak_col,
        "gene_column": gene_col,
    },
    "gene_activity": {
        "use_layer": use_layer,
        "normalize_total": bool(gene_activity_cfg.get("normalize_total", True)),
        "target_sum": float(gene_activity_cfg.get("target_sum", 10000.0)),
        "log1p": bool(gene_activity_cfg.get("log1p", True)),
        "min_peaks_per_gene": int(min_peaks_per_gene),
    },
    "mapping_summary": {
        "n_atac_peaks": n_atac_peaks,
        "n_annotation_rows": n_annotation_rows,
        "n_matched_peak_gene_edges": n_matched_peak_gene_edges,
        "n_matched_unique_peaks": n_matched_peaks,
        "fraction_atac_peaks_matched": float(n_matched_peaks / n_atac_peaks),
        "n_gene_activity_genes": n_genes,
    },
    "markers": {
        "n_requested_markers": len(marker_genes),
        "n_available_markers": len(available_markers),
        "available_markers": available_markers,
        "missing_markers": missing_markers,
    },
    "outputs": {
        "gene_activity_h5ad": str(gene_activity_out),
        "mapping_summary": str(mapping_summary_out),
        "marker_availability": str(marker_availability_out),
        "marker_mean_by_cluster": str(marker_mean_by_cluster_out),
        "marker_dotplot": str(marker_dotplot_out),
        "marker_umap": str(marker_umap_out),
    },
    "rules": {
        "outputs_project_local_only": True,
        "clustered_atac_object_not_overwritten": True,
        "rna_annotation_not_replaced": True,
        "gene_activity_used_for_atac_annotation_support": True,
    },
}

write_json(manifest, manifest_out)  # save manifest

print("\nSaved Step 04C manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 04C scATAC gene activity construction completed.")



print("\nExample ATAC peak IDs:")
print(list(atac.var_names[:5]))

print("\nExample annotation peak IDs:")
print(peak_annotation[peak_col].head().tolist())