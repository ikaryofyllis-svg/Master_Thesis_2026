"""Βήμα 03A — RNA marker validation ανά Leiden cluster"""

from pathlib import Path  # χειρισμός paths
import sys  # import path setup

import pandas as pd  # πίνακες για marker expression και availability
import scanpy as sc  # single-cell analysis και plots
import matplotlib.pyplot as plt  # save figures


# ============================================================
# 0. Ορισμός Refactored root και Python import path
# ============================================================

refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root

src_dir = refactor_dir / "src"  # shared package source folder

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # κάνει διαθέσιμο το pbmcgrn package


# ============================================================
# 1. Import helpers
# ============================================================

from pbmcgrn.config import load_yaml  # YAML loader
from pbmcgrn.io import write_csv, write_json  # save helpers
from pbmcgrn.paths import ProjectPaths  # project-aware paths


# ============================================================
# 2. Initialize project paths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # διαβάζει active project από registry
paths.ensure_output_dirs()  # δημιουργεί project folders αν λείπουν

print("\nActive project:")
print(paths.project_id)  # πρέπει να είναι pbmc_10k_v1

print("\nProject root:")
print(paths.project_root)  # πρέπει να είναι μέσα στο Refactored/projects/pbmc_10k_v1


# ============================================================
# 3. Load marker validation config
# ============================================================

config_path = paths.configs / "rna_marker_validation.yaml"  # project-local config
config = load_yaml(config_path)  # διαβάζει marker genes και outputs

cluster_key = config["cluster_key"]  # π.χ. leiden_0_6
marker_dict = config["marker_genes"]  # marker groups
outputs = config["outputs"]  # output paths

print("\nLoaded marker validation config:")
print(config_path)


# ============================================================
# 4. Helper functions
# ============================================================

def flatten_marker_dict(marker_dict):
    """
    Convert marker gene dictionary to one unique marker list.
    """

    markers = []  # εδώ μαζεύουμε όλα τα marker genes

    for genes in marker_dict.values():
        markers.extend(genes)  # προσθέτει markers από κάθε cell-type group

    markers = list(dict.fromkeys(markers))  # κρατά μοναδικά markers με αρχική σειρά

    return markers  # επιστρέφει flat marker list


def split_available_missing_markers(adata, marker_genes, use_raw=True):
    """
    Split marker genes into available and missing.
    """

    if use_raw and adata.raw is not None:
        var_names = adata.raw.var_names  # marker genes ψάχνονται στο raw full-gene matrix
    else:
        var_names = adata.var_names  # fallback στο τρέχον HVG matrix

    available = [gene for gene in marker_genes if gene in var_names]  # markers που υπάρχουν
    missing = [gene for gene in marker_genes if gene not in var_names]  # markers που λείπουν

    return available, missing  # επιστρέφει δύο λίστες


def compute_mean_marker_expression(adata, marker_genes, cluster_key, use_raw=True):
    """
    Compute mean marker expression per Leiden cluster.
    """

    available, missing = split_available_missing_markers(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )  # ελέγχει ποια markers υπάρχουν

    if len(available) == 0:
        raise ValueError("No marker genes were found in the AnnData object.")  # σταματά αν δεν βρεθεί τίποτα

    if use_raw and adata.raw is not None:
        expr_adata = adata.raw[:, available].to_adata()  # παίρνει normalized/log full-gene raw matrix
    else:
        expr_adata = adata[:, available].copy()  # παίρνει markers από current object

    expr_df = expr_adata.to_df()  # μετατρέπει expression matrix σε DataFrame
    expr_df[cluster_key] = adata.obs[cluster_key].astype(str).values  # προσθέτει Leiden cluster

    mean_expr = expr_df.groupby(cluster_key).mean()  # μέσο expression ανά cluster
    mean_expr = mean_expr.sort_index()  # σταθερή σειρά clusters

    return mean_expr, available, missing  # επιστρέφει table και marker lists


# ============================================================
# 5. Prepare output paths
# ============================================================

marker_dotplot_out = paths.prepare_output_path(
    outputs["marker_dotplot"],
    label="marker_dotplot",
)  # project-local marker dotplot

marker_umap_out = paths.prepare_output_path(
    outputs["marker_umap"],
    label="marker_umap",
)  # project-local marker UMAP

mean_expr_out = paths.prepare_output_path(
    outputs["marker_mean_expression_table"],
    label="marker_mean_expression_table",
)  # project-local mean marker expression table

marker_availability_out = paths.prepare_output_path(
    outputs["marker_availability_table"],
    label="marker_availability_table",
)  # project-local marker availability table

manifest_out = paths.prepare_output_path(
    outputs["manifest"],
    label="step03A_manifest",
)  # project-local manifest


# ============================================================
# 6. Load clustered RNA object
# ============================================================

rna_clustered_path = paths.processed_rna / "pbmc_rna_clustered.h5ad"  # Step 02 output

print("\nLoading clustered RNA object:")
print(rna_clustered_path)

if not rna_clustered_path.exists():
    raise FileNotFoundError(
        f"Clustered RNA object not found:\n{rna_clustered_path}\n"
        "Πρώτα πρέπει να τρέξει το Step 02."
    )  # σταματά αν δεν υπάρχει Step 02 output

rna = sc.read_h5ad(rna_clustered_path)  # φορτώνει clustered RNA object

print("\nRNA clustered object:")
print(rna)

if cluster_key not in rna.obs.columns:
    raise ValueError(f"Cluster key not found in rna.obs: {cluster_key}")  # ελέγχει ότι υπάρχει leiden_0_6


# ============================================================
# 7. Marker gene availability
# ============================================================

marker_genes = flatten_marker_dict(marker_dict)  # κάνει flat marker list

available_markers, missing_markers = split_available_missing_markers(
    adata=rna,
    marker_genes=marker_genes,
    use_raw=True,
)  # ελέγχει markers στο raw full-gene matrix αν υπάρχει

marker_availability = pd.DataFrame(
    {
        "marker_gene": marker_genes,
        "available": [gene in available_markers for gene in marker_genes],
    }
)  # table με available/missing markers

write_csv(marker_availability, marker_availability_out, index=False)  # σώζει marker availability

print("\nAvailable markers:")
print(available_markers)

print("\nMissing markers:")
print(missing_markers)

print("\nSaved marker availability table to:")
print(marker_availability_out)


# ============================================================
# 8. Marker dotplot by Leiden cluster
# ============================================================

print("\nCreating marker dotplot by Leiden cluster...")

sc.pl.dotplot(
    rna,
    marker_dict,
    groupby=cluster_key,
    use_raw=True if rna.raw is not None else False,
    standard_scale="var",
    show=False,
)  # dotplot markers ανά Leiden cluster

plt.savefig(
    marker_dotplot_out,
    dpi=300,
    bbox_inches="tight",
)  # αποθήκευση dotplot

plt.close()  # κλείνει figure

print("\nSaved marker dotplot to:")
print(marker_dotplot_out)


# ============================================================
# 9. Marker UMAPs
# ============================================================

print("\nCreating marker gene UMAPs...")

sc.pl.umap(
    rna,
    color=available_markers,
    use_raw=True if rna.raw is not None else False,
    frameon=False,
    show=False,
)  # UMAPs για κάθε marker gene

plt.savefig(
    marker_umap_out,
    dpi=300,
    bbox_inches="tight",
)  # αποθήκευση marker UMAP panel

plt.close()  # κλείνει figure

print("\nSaved marker UMAPs to:")
print(marker_umap_out)


# ============================================================
# 10. Mean marker expression by Leiden cluster
# ============================================================

print("\nComputing mean marker expression by Leiden cluster...")

mean_expr, available_markers, missing_markers = compute_mean_marker_expression(
    adata=rna,
    marker_genes=marker_genes,
    cluster_key=cluster_key,
    use_raw=True,
)  # μέσο marker expression ανά Leiden cluster

mean_expr.to_csv(mean_expr_out)  # αποθήκευση table

print("\nMean marker expression preview:")
print(mean_expr.round(2).head())

print("\nSaved mean marker expression table to:")
print(mean_expr_out)


# ============================================================
# 11. Save manifest
# ============================================================

manifest = {
    "step": "03A",
    "description": "RNA marker validation by Leiden cluster using marker dotplot, marker UMAPs and mean marker expression table",
    "project": paths.as_dict(),
    "config": str(config_path),
    "inputs": {
        "rna_clustered_h5ad": str(rna_clustered_path),
    },
    "outputs": {
        "marker_dotplot": str(marker_dotplot_out),
        "marker_umap": str(marker_umap_out),
        "marker_mean_expression_table": str(mean_expr_out),
        "marker_availability_table": str(marker_availability_out),
    },
    "summary": {
        "cluster_key": cluster_key,
        "n_cells": int(rna.n_obs),
        "n_genes_current_matrix": int(rna.n_vars),
        "has_raw": bool(rna.raw is not None),
        "n_marker_genes_requested": int(len(marker_genes)),
        "n_marker_genes_available": int(len(available_markers)),
        "n_marker_genes_missing": int(len(missing_markers)),
        "available_markers": available_markers,
        "missing_markers": missing_markers,
    },
}  # manifest για reproducibility

write_json(manifest, manifest_out)  # σώζει manifest

print("\nSaved manifest to:")
print(manifest_out)


# ============================================================
# 12. Final summary
# ============================================================

print("\nDone.")
print("Το Step 03A RNA marker validation ολοκληρώθηκε.")
print("Δεν προστέθηκαν final labels. Δημιουργήθηκαν μόνο marker evidence outputs.")