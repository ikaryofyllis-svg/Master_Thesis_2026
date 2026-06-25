"""Βήμα 03 — RNA marker validation, manual cell-type annotation και 5-type subset"""

from pathlib import Path
import sys

import scanpy as sc
import matplotlib.pyplot as plt


# ============================================================
# 0. Ορισμός Refactored project root και Python import path
# ============================================================

refactor_dir = Path(
    r"C:/Users/jkary/Desktop/bioinformatics/MSc thesis/Msc_Thesis_Project/Refactored"
)

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Imports από refactored package
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths

from pbmcgrn.preprocessing.annotation import (
    flatten_marker_dict,
    get_available_markers,
    compute_mean_marker_expression,
    apply_cluster_annotation,
    summarize_annotation,
    summarize_cell_type_counts,
    subset_to_scmultiomegrn_5types,
)


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("Refactored project root:")
print(paths.refactor_root)


# ============================================================
# 2. Load RNA annotation config
# ============================================================

config_path = paths.configs / "rna_annotation.yaml"
config = load_yaml(config_path)

cluster_key = config["cluster_key"]
marker_dict = config["marker_genes"]
detailed_mapping = config["cluster_to_cell_type_detailed"]
broad_mapping = config["cluster_to_cell_type_broad"]
scmultiomegrn_mapping = config["scmultiomegrn_5type_mapping"]
outputs = config["outputs"]

print("\nLoaded RNA annotation config:")
print(config_path)


# ============================================================
# 3. Helper function για paths από config
# ============================================================

def resolve_from_refactor(relative_or_absolute_path):
    """
    Resolve paths written inside config files.

    Relative paths are interpreted from Refactored/.
    """

    path = Path(relative_or_absolute_path)

    if path.is_absolute():
        return path

    return (paths.refactor_root / path).resolve()


# ============================================================
# 4. Load clustered RNA object
# ============================================================

rna_clustered_path = paths.processed_rna / "pbmc_rna_clustered.h5ad"

print("\nLoading clustered RNA object:")
print(rna_clustered_path)

if not rna_clustered_path.exists():
    raise FileNotFoundError(
        f"Clustered RNA object not found:\n{rna_clustered_path}\n"
        "Πρώτα πρέπει να τρέξει το Step 02 RNA preprocessing."
    )

rna = sc.read_h5ad(rna_clustered_path)

print("\nRNA clustered object:")
print(rna)

if cluster_key not in rna.obs.columns:
    raise ValueError(f"Cluster key not found in RNA object: {cluster_key}")


# ============================================================
# 5. Marker gene availability
# ============================================================

marker_genes = flatten_marker_dict(marker_dict)

available_markers, missing_markers = get_available_markers(
    adata=rna,
    marker_genes=marker_genes,
    use_raw=True,
)

print("\nAvailable marker genes:")
print(available_markers)

print("\nMissing marker genes:")
print(missing_markers)


# ============================================================
# 6. Marker dotplot by cluster
# ============================================================

marker_dotplot_out = resolve_from_refactor(outputs["marker_dotplot"])
marker_dotplot_out.parent.mkdir(parents=True, exist_ok=True)

sc.pl.dotplot(
    rna,
    marker_dict,
    groupby=cluster_key,
    use_raw=True if rna.raw is not None else False,
    standard_scale="var",
    show=False,
)

plt.savefig(
    marker_dotplot_out,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

print("\nSaved RNA marker dotplot to:")
print(marker_dotplot_out)


# ============================================================
# 7. Mean marker expression by cluster
# ============================================================

mean_expr, available_markers, missing_markers = compute_mean_marker_expression(
    adata=rna,
    marker_genes=marker_genes,
    cluster_key=cluster_key,
    use_raw=True,
)

mean_expr_out = resolve_from_refactor(outputs["marker_mean_expression_table"])
mean_expr_out.parent.mkdir(parents=True, exist_ok=True)

mean_expr.to_csv(mean_expr_out)

print("\nSaved mean marker expression table to:")
print(mean_expr_out)

print("\nMean marker expression preview:")
print(mean_expr.round(2).head())


# ============================================================
# 8. Apply manual cluster annotation
# ============================================================

print("\nApplying manual RNA cluster annotation...")

rna_annotated = apply_cluster_annotation(
    adata=rna,
    cluster_key=cluster_key,
    detailed_mapping=detailed_mapping,
    broad_mapping=broad_mapping,
)

annotation_summary = summarize_annotation(
    adata=rna_annotated,
    cluster_key=cluster_key,
)

cell_type_counts = summarize_cell_type_counts(rna_annotated)

print("\nAnnotation summary:")
print(annotation_summary)

print("\nCell type counts:")
print(cell_type_counts.head(20))


# ============================================================
# 9. Save annotated RNA object and annotation tables
# ============================================================

annotated_rna_out = resolve_from_refactor(outputs["annotated_rna_h5ad"])
annotated_rna_out.parent.mkdir(parents=True, exist_ok=True)

rna_annotated.write_h5ad(annotated_rna_out)

annotation_summary_out = resolve_from_refactor(
    outputs["cluster_annotation_summary"]
)

cell_type_counts_out = resolve_from_refactor(
    outputs["cell_type_counts"]
)

write_csv(annotation_summary, annotation_summary_out, index=False)
write_csv(cell_type_counts, cell_type_counts_out, index=False)

print("\nSaved annotated RNA object to:")
print(annotated_rna_out)

print("\nSaved annotation summary to:")
print(annotation_summary_out)

print("\nSaved cell type counts to:")
print(cell_type_counts_out)


# ============================================================
# 10. Save annotated UMAP plots
# ============================================================

umap_broad_out = resolve_from_refactor(outputs["annotated_umap_broad"])
umap_broad_out.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    rna_annotated,
    color="cell_type_broad",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    show=False,
)

plt.savefig(
    umap_broad_out,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

print("\nSaved broad cell-type UMAP to:")
print(umap_broad_out)


umap_detailed_out = resolve_from_refactor(outputs["annotated_umap_detailed"])
umap_detailed_out.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    rna_annotated,
    color="cell_type_detailed",
    legend_loc="right margin",
    frameon=False,
    show=False,
)

plt.savefig(
    umap_detailed_out,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

print("\nSaved detailed cell-type UMAP to:")
print(umap_detailed_out)


# ============================================================
# 11. Create 5-type RNA subset for downstream GRN pipeline
# ============================================================

print("\nCreating scMultiomeGRN-style 5-type RNA subset...")

rna_5, rna_5_counts = subset_to_scmultiomegrn_5types(
    adata=rna_annotated,
    mapping=scmultiomegrn_mapping,
)

rna_5_out = resolve_from_refactor(outputs["rna_5types_h5ad"])
rna_5_out.parent.mkdir(parents=True, exist_ok=True)

rna_5.write_h5ad(rna_5_out)

rna_5_counts_out = resolve_from_refactor(outputs["rna_5types_counts"])
write_csv(rna_5_counts, rna_5_counts_out, index=False)

print("\nRNA 5-type subset:")
print(rna_5)

print("\nRNA 5-type counts:")
print(rna_5_counts)

print("\nSaved RNA 5-type subset to:")
print(rna_5_out)

print("\nSaved RNA 5-type counts to:")
print(rna_5_counts_out)


# ============================================================
# 12. UMAP for 5-type subset
# ============================================================

rna_5_umap_out = resolve_from_refactor(outputs["rna_5types_umap"])
rna_5_umap_out.parent.mkdir(parents=True, exist_ok=True)

sc.pl.umap(
    rna_5,
    color="scmultiomegrn_cell_type",
    legend_loc="on data",
    legend_fontsize=6,
    legend_fontweight="normal",
    frameon=False,
    show=False,
)

plt.savefig(
    rna_5_umap_out,
    dpi=300,
    bbox_inches="tight",
)

plt.close()

print("\nSaved RNA 5-type UMAP to:")
print(rna_5_umap_out)


# ============================================================
# 13. Save manifest
# ============================================================

manifest = {
    "step": "03",
    "description": "RNA marker validation, manual cell-type annotation and 5-type RNA subset",
    "config": str(config_path),
    "inputs": {
        "rna_clustered_h5ad": str(rna_clustered_path),
    },
    "outputs": {
        "annotated_rna_h5ad": str(annotated_rna_out),
        "rna_5types_h5ad": str(rna_5_out),
        "marker_dotplot": str(marker_dotplot_out),
        "mean_marker_expression_table": str(mean_expr_out),
        "annotation_summary": str(annotation_summary_out),
        "cell_type_counts": str(cell_type_counts_out),
        "umap_broad": str(umap_broad_out),
        "umap_detailed": str(umap_detailed_out),
        "rna_5types_counts": str(rna_5_counts_out),
        "rna_5types_umap": str(rna_5_umap_out),
    },
    "summary": {
        "n_cells_annotated": int(rna_annotated.n_obs),
        "n_genes_annotated": int(rna_annotated.n_vars),
        "n_cells_5type_subset": int(rna_5.n_obs),
        "n_genes_5type_subset": int(rna_5.n_vars),
        "available_markers": available_markers,
        "missing_markers": missing_markers,
    },
}

manifest_out = paths.manifests / "step03_annotate_rna_manifest.json"
write_json(manifest, manifest_out)

print("\nSaved manifest to:")
print(manifest_out)


# ============================================================
# 14. Τελικό summary
# ============================================================

print("\nDone.")
print("Το Step 03 RNA annotation ολοκληρώθηκε.")
print("Έχουμε annotated RNA object και 5-type RNA subset για downstream integration/GRN.")