
# ============================================================
# ATAC clustering stability check using atac_lsi
# ============================================================

import scanpy as sc  # single-cell analysis toolkit
import pandas as pd  # tables / summaries
from pathlib import Path  # safe path handling

# ------------------------------------------------------------
# 1. Load clustered ATAC object as atac_lsi
# ------------------------------------------------------------

cwd = Path.cwd().resolve()  # current notebook/script folder

if cwd.name == "Refactored":  # if notebook is already inside Refactored/
    refactor_dir = cwd  # use current folder as project root
else:
    refactor_dir = cwd / "Refactored"  # otherwise expect Refactored/ inside cwd

active_project = refactor_dir / "projects" / "pbmc_10k_v1"  # active project folder

atac_h5ad = active_project / "data" / "processed" / "atac" / "pbmc_atac_clustered.h5ad"  # clustered ATAC object

print("ATAC file exists?", atac_h5ad.exists())  # check file existence
print("ATAC path:", atac_h5ad)  # show full path

# -*- coding: utf-8 -*-
"""
Step 04B2 — exploratory scATAC clustering stability check.

Purpose:
- Load the official Step 04B clustered ATAC object.
- Re-run Leiden clustering across several n_neighbors and resolution values.
- Save a stability summary table and UMAP comparison figure.
- Do not overwrite the official clustered ATAC object.
"""

from pathlib import Path
import sys
from datetime import datetime

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root
except NameError:
    refactor_dir = Path.cwd().resolve()  # interactive fallback

src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # allow local pbmcgrn imports

from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


def save_show_close(path: Path, dpi: int = 300, show: bool = True) -> None:
    """Save current matplotlib figure and optionally show it."""

    plt.savefig(path, dpi=dpi, bbox_inches="tight")  # save figure

    if show:
        plt.show()  # display figure

    plt.close()  # close figure


paths = ProjectPaths(refactor_root=refactor_dir)  # project-aware paths
paths.ensure_output_dirs()  # create standard output dirs

print("\nActive project:")
print(paths.project_id)

print("\nProject root:")
print(paths.project_root)


# ============================================================
# 1. Resolve inputs / outputs
# ============================================================

clustered_in = paths.project_root / "data/processed/atac/pbmc_atac_clustered.h5ad"  # official 04B output

table_out = paths.project_root / "reports/tables/atac/atac_clustering_stability_exploratory.csv"  # stability summary

fig_out = paths.project_root / "reports/figures/atac/clustering_stability/atac_umap_clustering_stability_candidates.png"  # UMAP comparison

manifest_out = paths.project_root / "reports/manifests/step04B2_scatac_clustering_stability_manifest.json"  # manifest

table_out.parent.mkdir(parents=True, exist_ok=True)  # create table folder
fig_out.parent.mkdir(parents=True, exist_ok=True)  # create figure folder
manifest_out.parent.mkdir(parents=True, exist_ok=True)  # create manifest folder

if not clustered_in.exists():
    raise FileNotFoundError(f"Clustered ATAC object not found:\n{clustered_in}")


# ============================================================
# 2. Load ATAC object and run sanity checks
# ============================================================

atac_lsi = sc.read_h5ad(clustered_in)  # load official clustered ATAC object

print("\nLoaded ATAC object:")
print(atac_lsi)

print("\nobsm keys:")
print(list(atac_lsi.obsm.keys()))

print("\nobs columns:")
print(list(atac_lsi.obs.columns))

if "X_lsi_2_50" not in atac_lsi.obsm.keys():
    raise KeyError("X_lsi_2_50 not found in atac_lsi.obsm. Check Step 04B output.")

if "leiden_0_6" not in atac_lsi.obs.columns:
    raise KeyError("leiden_0_6 not found in atac_lsi.obs. Check Step 04B output.")

print("\nCurrent official ATAC clusters:")
print(atac_lsi.obs["leiden_0_6"].value_counts().sort_index())


# ============================================================
# 3. Stability screen: neighbors × Leiden resolution
# ============================================================

results = []  # collect one row per clustering setting

for n_neighbors in [10, 15, 20, 30]:  # test local-vs-global graph structure

    neighbors_key = f"neighbors_nn{n_neighbors}"  # graph key for this n_neighbors

    sc.pp.neighbors(
        atac_lsi,
        n_neighbors=n_neighbors,
        use_rep="X_lsi_2_50",
        metric="cosine",
        key_added=neighbors_key,
    )  # build neighbor graph from ATAC LSI space

    for resolution in [0.3, 0.4, 0.6, 0.8, 1.0]:  # test coarser/finer Leiden settings

        cluster_key = f"leiden_nn{n_neighbors}_res{str(resolution).replace('.', '_')}"  # safe column name

        sc.tl.leiden(
            atac_lsi,
            resolution=resolution,
            key_added=cluster_key,
            neighbors_key=neighbors_key,
        )  # run Leiden on this graph

        counts = atac_lsi.obs[cluster_key].value_counts().sort_index()  # cluster sizes

        results.append(
            {
                "n_neighbors": int(n_neighbors),
                "resolution": float(resolution),
                "cluster_key": cluster_key,
                "n_clusters": int(counts.shape[0]),
                "min_cluster_size": int(counts.min()),
                "max_cluster_size": int(counts.max()),
                "median_cluster_size": float(counts.median()),
            }
        )  # save summary row

        print(cluster_key, counts.to_dict())  # print cluster sizes


stability_summary = pd.DataFrame(results)  # convert diagnostics to table

stability_summary = stability_summary.sort_values(
    ["n_neighbors", "resolution"]
).reset_index(drop=True)  # clean ordering

write_csv(stability_summary, table_out, index=False)  # save audit table

print("\nSaved stability summary to:")
print(table_out)


# ============================================================
# 4. UMAP comparison for useful candidate clusterings
# ============================================================

candidate_keys = [
    "leiden_0_6",
    "leiden_nn15_res0_4",
    "leiden_nn20_res0_4",
    "leiden_nn20_res0_6",
]  # selected clusterings for visual comparison

missing_keys = [key for key in candidate_keys if key not in atac_lsi.obs.columns]  # check required columns

if missing_keys:
    raise KeyError(f"Missing candidate cluster keys: {missing_keys}")

sc.pl.umap(
    atac_lsi,
    color=candidate_keys,
    legend_loc="on data",
    ncols=2,
    wspace=0.4,
    show=False,
)  # compare official vs coarser/smoother alternatives

save_show_close(fig_out, dpi=300, show=True)  # save comparison figure

print("\nSaved UMAP comparison figure to:")
print(fig_out)


# ============================================================
# 5. Manifest
# ============================================================

manifest = {
    "step": "04B2",
    "description": "Exploratory scATAC clustering stability screen across neighbors and Leiden resolutions",
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "inputs": {
        "clustered_atac_h5ad": str(clustered_in),
    },
    "parameters": {
        "n_neighbors_values": [10, 15, 20, 30],
        "resolution_values": [0.3, 0.4, 0.6, 0.8, 1.0],
        "use_rep": "X_lsi_2_50",
        "metric": "cosine",
    },
    "outputs": {
        "stability_summary": str(table_out),
        "umap_comparison": str(fig_out),
    },
    "interpretation": {
        "official_clustering_retained": True,
        "official_cluster_key": "leiden_0_6",
        "purpose": "diagnostic only; does not overwrite official Step 04B output",
    },
}

write_json(manifest, manifest_out)  # save manifest

print("\nSaved Step 04B2 manifest to:")
print(manifest_out)

print("\nDone.")
print("Step 04B2 scATAC clustering stability check completed.")
