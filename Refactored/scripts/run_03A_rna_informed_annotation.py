"""
Step 03A — RNA informed marker-based annotation evidence.

Purpose:
- Load the final selected clustered RNA object from Step 02B.
- Read the selected final Leiden cluster key dynamically from preprocessing_rna.yaml.
- Generate marker-based evidence tables and plots for manual/informed annotation.
- Do NOT create final annotation labels.
- Do NOT write final annotated h5ad.
- Do NOT rerun QC, PCA, neighbors, UMAP, or Leiden.
- Write all outputs inside the active project only.

Source of truth for clustering:
configs/preprocessing_rna.yaml
  -> dimensionality_reduction
     -> leiden_key
     -> n_pcs
     -> n_neighbors
     -> leiden_resolution
"""

from __future__ import annotations  # cleaner type hints

from pathlib import Path  # robust filesystem paths
import sys  # allows adding Refactored/src to Python import path
from typing import Any  # flexible typing for YAML dictionaries
from datetime import datetime  # manifest timestamp
import re  # safe filename sanitization

import pandas as pd  # evidence tables
import scanpy as sc  # single-cell analysis and plotting
import matplotlib.pyplot as plt  # save Scanpy figures


# ============================================================
# 0. Resolve Refactored root and Python import path
# ============================================================

try:
    refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root when script is in Refactored/scripts
except NameError:
    refactor_dir = Path.cwd().resolve()  # fallback for interactive execution

src_dir = refactor_dir / "src"  # shared reusable source code folder

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))  # makes pbmcgrn importable


# ============================================================
# 1. Import project helpers
# ============================================================

from pbmcgrn.config import load_yaml  # YAML loader
from pbmcgrn.io import write_csv, write_json  # reusable output writers
from pbmcgrn.paths import ProjectPaths  # project-aware path manager


# ============================================================
# 2. Generic validation and path helpers
# ============================================================

def assert_path_inside_project(path: Path, project_root: Path, label: str) -> None:
    """Fail if a path is outside the active project root."""

    resolved_path = Path(path).resolve()  # normalize path
    resolved_project_root = Path(project_root).resolve()  # normalize active project root

    try:
        resolved_path.relative_to(resolved_project_root)  # succeeds only if path is inside project
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project workspace.\n"
            f"Path: {resolved_path}\n"
            f"Active project root: {resolved_project_root}"
        ) from exc  # fail safely


def require_config_key(config: dict[str, Any], key: str, config_name: str) -> Any:
    """Return required config value or fail clearly."""

    if key not in config:
        raise KeyError(
            f"Missing required key '{key}' in {config_name}."
        )  # avoids silent config drift

    return config[key]  # required value


def ensure_parent_dir(path: Path) -> None:
    """Create parent directory for an output file."""

    Path(path).parent.mkdir(parents=True, exist_ok=True)  # ensures save location exists
    
def save_and_optionally_show_figure(
    output_path: Path,
    dpi: int,
    display_figures: bool,
) -> None:
    """Save current matplotlib figure and optionally display it."""

    plt.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )  # save figure to project-local output path

    if display_figures:
        plt.show()  # display figure in Spyder/Jupyter/interactive backend
    else:
        plt.close()  # close figure only when not displaying

    if display_figures:
        plt.close()  # close after display to avoid memory buildup

def safe_filename_token(text: str) -> str:
    """Convert a config key into a safe filename token."""

    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text))  # replace unsafe filename chars

    return safe.strip("_")  # remove leading/trailing underscores


def sort_cluster_table(df: pd.DataFrame, cluster_col: str = "cluster") -> pd.DataFrame:
    """Sort clusters numerically when possible, otherwise lexicographically."""

    sorted_df = df.copy()  # avoid mutating original table

    sorted_df["_cluster_sort"] = pd.to_numeric(
        sorted_df[cluster_col],
        errors="coerce",
    )  # numeric sort helper for Leiden clusters like 0, 1, 2, ...

    sorted_df = sorted_df.sort_values(
        by=["_cluster_sort", cluster_col],
        na_position="last",
    )  # stable cluster ordering

    sorted_df = sorted_df.drop(columns=["_cluster_sort"])  # remove helper column

    return sorted_df  # sorted table


# ============================================================
# 3. Marker helper functions
# ============================================================

def flatten_marker_dict(marker_dict: dict[str, list[str]]) -> list[str]:
    """Flatten marker dictionary into one unique ordered marker list."""

    markers: list[str] = []  # stores all requested marker genes

    for genes in marker_dict.values():
        markers.extend(genes)  # append marker genes from each biological group

    markers = list(dict.fromkeys(markers))  # remove duplicates while preserving order

    return markers  # unique marker list


def get_marker_var_names(adata: sc.AnnData, use_raw: bool) -> pd.Index:
    """Return gene names used for marker availability checking."""

    if use_raw and adata.raw is not None:
        return adata.raw.var_names  # preferred: full raw/log-normalized gene matrix

    return adata.var_names  # fallback: current matrix, often HVG-only


def split_available_missing_markers(
    adata: sc.AnnData,
    marker_genes: list[str],
    use_raw: bool,
) -> tuple[list[str], list[str]]:
    """Split requested marker genes into available and missing."""

    var_names = get_marker_var_names(
        adata=adata,
        use_raw=use_raw,
    )  # gene universe to check against

    available = [
        gene for gene in marker_genes
        if gene in var_names
    ]  # markers found in AnnData

    missing = [
        gene for gene in marker_genes
        if gene not in var_names
    ]  # markers missing from AnnData

    return available, missing  # marker availability split


def filter_marker_dict_to_available(
    marker_dict: dict[str, list[str]],
    available_markers: list[str],
) -> dict[str, list[str]]:
    """Keep only available markers in each marker group."""

    available_set = set(available_markers)  # fast membership checks

    filtered = {
        group: [
            gene for gene in genes
            if gene in available_set
        ]
        for group, genes in marker_dict.items()
    }  # preserve biological marker groups

    filtered = {
        group: genes
        for group, genes in filtered.items()
        if len(genes) > 0
    }  # remove empty marker groups

    return filtered  # Scanpy-safe marker dictionary


def build_marker_availability_table(
    marker_dict: dict[str, list[str]],
    available_markers: list[str],
) -> pd.DataFrame:
    """Create long marker availability table."""

    available_set = set(available_markers)  # fast membership checks

    rows = []  # table rows

    for marker_group, genes in marker_dict.items():
        for gene in genes:
            rows.append(
                {
                    "marker_group": marker_group,  # biological marker group
                    "marker_gene": gene,  # requested marker gene
                    "available": gene in available_set,  # found or missing
                }
            )  # one row per marker occurrence

    availability = pd.DataFrame(rows)  # long-format marker audit table

    return availability  # output table


def get_expression_dataframe(
    adata: sc.AnnData,
    marker_genes: list[str],
    use_raw: bool,
) -> pd.DataFrame:
    """Extract marker expression matrix as cells x marker genes DataFrame."""

    if len(marker_genes) == 0:
        raise ValueError(
            "No available marker genes were provided."
        )  # cannot compute marker summaries

    if use_raw and adata.raw is not None:
        expr_adata = adata.raw[:, marker_genes].to_adata()  # use full raw/log matrix
    else:
        expr_adata = adata[:, marker_genes].copy()  # use current matrix fallback

    expr_df = expr_adata.to_df()  # convert to pandas DataFrame

    return expr_df  # cell-level expression table


def compute_mean_marker_expression(
    adata: sc.AnnData,
    marker_genes: list[str],
    cluster_key: str,
    use_raw: bool,
) -> pd.DataFrame:
    """Compute mean marker expression per cluster."""

    expr_df = get_expression_dataframe(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )  # cells x markers

    expr_df[cluster_key] = adata.obs[cluster_key].astype(str).values  # attach cluster labels

    mean_expr = (
        expr_df
        .groupby(cluster_key)
        .mean(numeric_only=True)
        .reset_index()
        .rename(columns={cluster_key: "cluster"})
    )  # cluster x mean marker expression

    mean_expr = sort_cluster_table(
        mean_expr,
        cluster_col="cluster",
    )  # natural cluster sorting

    return mean_expr  # output table


def compute_pct_marker_expressing(
    adata: sc.AnnData,
    marker_genes: list[str],
    cluster_key: str,
    use_raw: bool,
) -> pd.DataFrame:
    """Compute percent of cells expressing each marker per cluster."""

    expr_df = get_expression_dataframe(
        adata=adata,
        marker_genes=marker_genes,
        use_raw=use_raw,
    )  # cells x markers

    binary_expr = (expr_df > 0).astype(float)  # detected expression indicator

    binary_expr[cluster_key] = adata.obs[cluster_key].astype(str).values  # attach cluster labels

    pct_expr = (
        binary_expr
        .groupby(cluster_key)
        .mean(numeric_only=True)
        .mul(100.0)
        .reset_index()
        .rename(columns={cluster_key: "cluster"})
    )  # cluster x percent expressing

    pct_expr = sort_cluster_table(
        pct_expr,
        cluster_col="cluster",
    )  # natural cluster sorting

    return pct_expr  # output table


def compute_cluster_sizes(
    adata: sc.AnnData,
    cluster_key: str,
) -> pd.DataFrame:
    """Compute cluster sizes and cell fractions."""

    cluster_sizes = (
        adata.obs[cluster_key]
        .astype(str)
        .value_counts()
        .rename_axis("cluster")
        .reset_index(name="n_cells")
    )  # one row per cluster

    cluster_sizes["fraction_cells"] = (
        cluster_sizes["n_cells"] / adata.n_obs
    )  # cluster fraction among all cells

    cluster_sizes = sort_cluster_table(
        cluster_sizes,
        cluster_col="cluster",
    )  # natural cluster sorting

    return cluster_sizes  # output table


def compute_marker_group_scores(
    mean_expr: pd.DataFrame,
    marker_dict_available: dict[str, list[str]],
) -> pd.DataFrame:
    """Compute simple average marker-group signal per cluster."""

    score_table = mean_expr[["cluster"]].copy()  # preserve cluster column

    for marker_group, genes in marker_dict_available.items():
        genes_in_table = [
            gene for gene in genes
            if gene in mean_expr.columns
        ]  # only genes present in mean expression table

        if len(genes_in_table) == 0:
            continue  # skip unavailable marker groups

        score_table[marker_group] = mean_expr[genes_in_table].mean(axis=1)  # group-level average signal

    return score_table  # cluster x marker group score table


def build_annotation_template(
    cluster_sizes: pd.DataFrame,
    marker_group_scores: pd.DataFrame,
) -> pd.DataFrame:
    """Build manual annotation worksheet without final labels."""

    score_columns = [
        col for col in marker_group_scores.columns
        if col != "cluster"
    ]  # marker group score columns

    if len(score_columns) == 0:
        top_signal = pd.DataFrame(
            {
                "cluster": cluster_sizes["cluster"].astype(str),
                "top_marker_signal": "",
            }
        )  # fallback if no marker scores exist
    else:
        score_copy = marker_group_scores.copy()  # avoid mutating original score table

        score_copy["top_marker_signal"] = score_copy[score_columns].idxmax(
            axis=1
        )  # strongest marker group by mean expression

        top_signal = score_copy[
            ["cluster", "top_marker_signal"]
        ]  # compact helper table

    template = cluster_sizes.merge(
        top_signal,
        on="cluster",
        how="left",
    )  # attach strongest marker group to each cluster

    template["candidate_cell_type"] = ""  # manual-fill candidate label, not final
    template["confidence"] = ""  # manual-fill confidence
    template["notes"] = ""  # manual notes
    template["needs_subclustering"] = ""  # manual flag for mixed/broad clusters

    template = template[
        [
            "cluster",
            "n_cells",
            "fraction_cells",
            "top_marker_signal",
            "candidate_cell_type",
            "confidence",
            "notes",
            "needs_subclustering",
        ]
    ]  # thesis-ready worksheet column order

    return template  # annotation evidence worksheet


def save_marker_umap_batches(
    adata: sc.AnnData,
    available_markers: list[str],
    marker_umaps_dir: Path,
    batch_size: int,
    use_raw: bool,
    dpi: int,
    cluster_key: str,
    display_figures: bool,
) -> list[str]:
    """Save marker UMAP panels in batches and optionally display them."""

    if "X_umap" not in adata.obsm:
        raise KeyError(
            "UMAP coordinates were not found in adata.obsm['X_umap'].\n"
            "Step 03A must not recompute UMAP. Rerun Step 02B first."
        )  # no UMAP evidence without Step 02B UMAP

    marker_umaps_dir.mkdir(parents=True, exist_ok=True)  # create output folder

    saved_files: list[str] = []  # list of saved files

    for batch_start in range(0, len(available_markers), batch_size):
        marker_batch = available_markers[
            batch_start:batch_start + batch_size
        ]  # marker subset for this panel

        batch_id = (batch_start // batch_size) + 1  # one-based batch index

        batch_out = marker_umaps_dir / (
            f"rna_marker_umaps_{cluster_key}_batch_{batch_id:02d}.png"
        )  # output figure path

        sc.pl.umap(
            adata,
            color=marker_batch,
            use_raw=use_raw,
            frameon=False,
            show=False,
        )  # create marker expression UMAP panel without immediate display

        save_and_optionally_show_figure(
            output_path=batch_out,
            dpi=dpi,
            display_figures=display_figures,
        )  # save and optionally show marker UMAP batch

        saved_files.append(str(batch_out))  # track output path

    return saved_files  # saved UMAP files

# ============================================================
# 4. Initialize project-aware paths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # active project path manager

paths.ensure_output_dirs()  # create standard project output folders

print("\nRefactored project root:")
print(paths.refactor_root)  # expected Refactored root

print("\nActive project id:")
print(paths.project_id)  # expected active project id

print("\nActive project root:")
print(paths.project_root)  # expected active project directory


# ============================================================
# 5. Load Step 03A marker config
# ============================================================

marker_config_path = paths.configs / "rna_marker_validation.yaml"  # marker evidence config

if not marker_config_path.exists():
    raise FileNotFoundError(
        f"RNA marker validation config not found:\n{marker_config_path}"
    )  # cannot run Step 03A without marker config

assert_path_inside_project(
    marker_config_path,
    paths.project_root,
    label="rna_marker_validation.yaml",
)  # config must be active-project-local

marker_config = load_yaml(marker_config_path)  # load Step 03A config

input_cfg = require_config_key(
    marker_config,
    "input",
    "rna_marker_validation.yaml",
)  # input section

marker_dict = require_config_key(
    marker_config,
    "marker_genes",
    "rna_marker_validation.yaml",
)  # canonical PBMC marker panel

outputs = require_config_key(
    marker_config,
    "outputs",
    "rna_marker_validation.yaml",
)  # output paths

plotting_cfg = marker_config.get(
    "plotting",
    {},
)  # plotting settings with defaults

use_raw_requested = bool(
    plotting_cfg.get("use_raw", True)
)  # prefer raw full-gene matrix if available

dotplot_standard_scale = plotting_cfg.get(
    "dotplot_standard_scale",
    "var",
)  # Scanpy marker plot scaling

umap_marker_batch_size = int(
    plotting_cfg.get("umap_marker_batch_size", 12)
)  # number of marker genes per UMAP panel

dpi = int(
    plotting_cfg.get("dpi", 300)
)  # figure resolution

display_figures = bool(
    plotting_cfg.get("display_figures", True)
)  # show figures interactively as well as saving them
print("\nLoaded Step 03A marker config:")
print(marker_config_path)


# ============================================================
# 6. Load Step 02B preprocessing config as clustering source of truth
# ============================================================

preprocessing_config_path = paths.configs / "preprocessing_rna.yaml"  # selected preprocessing config

if not preprocessing_config_path.exists():
    raise FileNotFoundError(
        f"RNA preprocessing config not found:\n{preprocessing_config_path}\n"
        "Step 03A needs this file to read the selected final Leiden key from Step 02B."
    )  # cannot infer selected cluster key

assert_path_inside_project(
    preprocessing_config_path,
    paths.project_root,
    label="preprocessing_rna.yaml",
)  # preprocessing config must be project-local

preprocessing_config = load_yaml(
    preprocessing_config_path
)  # load Step 02B selected preprocessing config

dim_cfg = require_config_key(
    preprocessing_config,
    "dimensionality_reduction",
    "preprocessing_rna.yaml",
)  # selected PCA/neighbors/Leiden config

required_dim_keys = [
    "n_pcs",
    "n_neighbors",
    "leiden_resolution",
    "leiden_key",
]  # keys required for Step 02B -> Step 03A handoff

missing_dim_keys = [
    key for key in required_dim_keys
    if key not in dim_cfg
]  # detect incomplete selected clustering config

if missing_dim_keys:
    raise KeyError(
        f"Missing required keys in preprocessing_rna.yaml::dimensionality_reduction: "
        f"{missing_dim_keys}"
    )  # fail clearly

cluster_key = str(
    dim_cfg["leiden_key"]
)  # dynamic selected cluster key from Step 02B

cluster_key_file_token = safe_filename_token(
    cluster_key
)  # safe token for dynamic output filenames

selected_n_pcs = int(
    dim_cfg["n_pcs"]
)  # final selected PCs

selected_n_neighbors = int(
    dim_cfg["n_neighbors"]
)  # final selected neighbors

selected_leiden_resolution = float(
    dim_cfg["leiden_resolution"]
)  # final selected Leiden resolution

print("\nStep 03A selected clustering source:")
print(f"source config: {preprocessing_config_path}")
print(f"cluster_key: {cluster_key}")
print(f"n_pcs: {selected_n_pcs}")
print(f"n_neighbors: {selected_n_neighbors}")
print(f"leiden_resolution: {selected_leiden_resolution}")


# ============================================================
# 7. Resolve project-local input path
# ============================================================

rna_clustered_path = paths.prepare_output_path(
    input_cfg["rna_clustered_h5ad"],
    label="rna_clustered_h5ad_input",
)  # project-local clustered RNA h5ad from Step 02B

assert_path_inside_project(
    rna_clustered_path,
    paths.project_root,
    label="rna_clustered_h5ad_input",
)  # input must be inside active project

if not rna_clustered_path.exists():
    raise FileNotFoundError(
        f"Clustered RNA object not found:\n{rna_clustered_path}\n\n"
        "Run Step 02B first and confirm that pbmc_rna_clustered.h5ad exists."
    )  # Step 03A depends on selected Step 02B object


# ============================================================
# 8. Resolve project-local output paths
# ============================================================

marker_availability_out = paths.prepare_output_path(
    outputs.get(
        "marker_availability_table",
        "reports/tables/rna/rna_marker_gene_availability.csv",
    ),
    label="marker_availability_table",
)  # marker availability audit table

mean_expr_out = paths.prepare_output_path(
    outputs.get(
        "marker_mean_expression_table",
        "reports/tables/rna/rna_cluster_marker_mean_expression.csv",
    ),
    label="marker_mean_expression_table",
)  # mean marker expression table

pct_expr_out = paths.prepare_output_path(
    outputs.get(
        "marker_pct_expressing_table",
        "reports/tables/rna/rna_cluster_marker_pct_expressing.csv",
    ),
    label="marker_pct_expressing_table",
)  # percent expressing table

marker_group_scores_out = paths.prepare_output_path(
    outputs.get(
        "marker_group_scores_table",
        "reports/tables/rna/rna_cluster_marker_group_scores.csv",
    ),
    label="marker_group_scores_table",
)  # marker-group score table

annotation_template_out = paths.prepare_output_path(
    outputs.get(
        "annotation_template_table",
        "reports/tables/rna/rna_informed_annotation_template.csv",
    ),
    label="annotation_template_table",
)  # manual annotation worksheet

cluster_sizes_out = paths.prepare_output_path(
    outputs.get(
        "cluster_sizes_table",
        "reports/tables/rna/rna_step03A_cluster_sizes.csv",
    ),
    label="cluster_sizes_table",
)  # cluster size audit table

marker_dotplot_dir = paths.prepare_output_path(
    outputs.get(
        "marker_dotplot_dir",
        "reports/figures/rna/annotation",
    ),
    label="marker_dotplot_dir",
)  # annotation figure directory

marker_umaps_base_dir = paths.prepare_output_path(
    outputs.get(
        "marker_umaps_dir",
        "reports/figures/rna/annotation/marker_umaps",
    ),
    label="marker_umaps_dir",
)  # base directory for marker UMAPs

marker_umaps_dir = marker_umaps_base_dir / cluster_key_file_token  # separate UMAP folder per cluster key

marker_dotplot_out = marker_dotplot_dir / (
    f"rna_marker_dotplot_by_{cluster_key_file_token}.png"
)  # dynamic dotplot output

marker_matrixplot_out = marker_dotplot_dir / (
    f"rna_marker_matrixplot_by_{cluster_key_file_token}.png"
)  # dynamic matrixplot output

manifest_out = paths.prepare_output_path(
    outputs.get(
        "manifest",
        "reports/manifests/step03A_rna_informed_annotation_manifest.json",
    ),
    label="step03A_manifest",
)  # Step 03A manifest

all_project_paths = {
    "rna_clustered_h5ad_input": rna_clustered_path,
    "marker_availability_table": marker_availability_out,
    "marker_mean_expression_table": mean_expr_out,
    "marker_pct_expressing_table": pct_expr_out,
    "marker_group_scores_table": marker_group_scores_out,
    "annotation_template_table": annotation_template_out,
    "cluster_sizes_table": cluster_sizes_out,
    "marker_dotplot_dir": marker_dotplot_dir,
    "marker_dotplot": marker_dotplot_out,
    "marker_matrixplot": marker_matrixplot_out,
    "marker_umaps_base_dir": marker_umaps_base_dir,
    "marker_umaps_dir": marker_umaps_dir,
    "step03A_manifest": manifest_out,
}  # every important path to validate

for label, path in all_project_paths.items():
    assert_path_inside_project(
        path,
        paths.project_root,
        label=label,
    )  # strict project-local validation

for output_file in [
    marker_availability_out,
    mean_expr_out,
    pct_expr_out,
    marker_group_scores_out,
    annotation_template_out,
    cluster_sizes_out,
    marker_dotplot_out,
    marker_matrixplot_out,
    manifest_out,
]:
    ensure_parent_dir(output_file)  # create output parent directories

marker_dotplot_dir.mkdir(
    parents=True,
    exist_ok=True,
)  # create figure output directory

marker_umaps_dir.mkdir(
    parents=True,
    exist_ok=True,
)  # create cluster-key-specific UMAP directory

print("\nResolved Step 03A input/output paths:")
for label, path in all_project_paths.items():
    print(f"{label}: {path}")  # terminal audit


# ============================================================
# 9. Load clustered RNA object
# ============================================================

print("\nLoading clustered RNA object:")
print(rna_clustered_path)

rna = sc.read_h5ad(
    rna_clustered_path
)  # load final selected clustered RNA object

print("\nRNA clustered object:")
print(rna)  # AnnData summary


# ============================================================
# 10. Validate selected cluster key and UMAP exist
# ============================================================

if cluster_key not in rna.obs.columns:
    available_obs_columns = list(rna.obs.columns)  # useful debugging output

    raise ValueError(
        f"The selected cluster key from preprocessing_rna.yaml was not found in rna.obs.\n\n"
        f"Selected cluster_key: {cluster_key}\n"
        f"Available rna.obs columns:\n{available_obs_columns}\n\n"
        "This usually means preprocessing_rna.yaml was updated, but Step 02B was not rerun, "
        "or pbmc_rna_clustered.h5ad is not the selected final object."
    )  # fail rather than using stale clusters

if "X_umap" not in rna.obsm:
    raise KeyError(
        "UMAP coordinates not found in rna.obsm['X_umap'].\n"
        "Step 03A must not recompute UMAP. Rerun Step 02B first."
    )  # Step 03A is evidence-only, not preprocessing

cluster_counts_preview = (
    rna.obs[cluster_key]
    .astype(str)
    .value_counts()
    .rename_axis("cluster")
    .reset_index(name="n_cells")
)  # cluster count preview

cluster_counts_preview = sort_cluster_table(
    cluster_counts_preview,
    cluster_col="cluster",
)  # natural cluster sorting

print("\nSelected cluster counts:")
print(cluster_counts_preview)  # terminal sanity check


# ============================================================
# 11. Marker availability
# ============================================================

use_raw = bool(
    use_raw_requested and rna.raw is not None
)  # actual raw usage

if use_raw_requested and rna.raw is None:
    print(
        "\nWARNING: plotting.use_raw=true, but rna.raw is missing. "
        "Falling back to adata.var_names/current matrix."
    )  # warning, not fatal

marker_genes = flatten_marker_dict(
    marker_dict
)  # unique requested markers for tables

available_markers, missing_markers = split_available_missing_markers(
    adata=rna,
    marker_genes=marker_genes,
    use_raw=use_raw,
)  # marker availability check

if len(available_markers) == 0:
    raise ValueError(
        "No requested marker genes were found in the RNA object.\n"
        "Check marker symbols, gene naming, raw matrix availability, and Step 02B output."
    )  # cannot continue without marker evidence

marker_dict_available = filter_marker_dict_to_available(
    marker_dict=marker_dict,
    available_markers=available_markers,
)  # marker dictionary filtered for Scanpy plots

marker_availability = build_marker_availability_table(
    marker_dict=marker_dict,
    available_markers=available_markers,
)  # marker availability audit table

write_csv(
    marker_availability,
    marker_availability_out,
    index=False,
)  # save marker availability table

print("\nMarker availability summary:")
print(f"Requested unique markers: {len(marker_genes)}")
print(f"Available unique markers: {len(available_markers)}")
print(f"Missing unique markers: {len(missing_markers)}")

print("\nMissing markers:")
print(missing_markers)

print("\nSaved marker availability table to:")
print(marker_availability_out)


# ============================================================
# 12. Cluster size table
# ============================================================

cluster_sizes = compute_cluster_sizes(
    adata=rna,
    cluster_key=cluster_key,
)  # cluster size evidence table

write_csv(
    cluster_sizes,
    cluster_sizes_out,
    index=False,
)  # save cluster sizes

print("\nSaved cluster sizes table to:")
print(cluster_sizes_out)


# ============================================================
# 13. Marker dotplot
# ============================================================

print("\nCreating marker dotplot by selected cluster key...")

sc.pl.dotplot(
    rna,
    marker_dict_available,
    groupby=cluster_key,
    use_raw=use_raw,
    standard_scale=dotplot_standard_scale,
    show=False,
)  # create dotplot without immediate display so we can save first

save_and_optionally_show_figure(
    output_path=marker_dotplot_out,
    dpi=dpi,
    display_figures=display_figures,
)  # save and optionally show dotplot

plt.close()  # release matplotlib figure

print("\nSaved marker dotplot to:")
print(marker_dotplot_out)


# ============================================================
# 14. Marker matrixplot
# ============================================================

print("\nCreating marker matrixplot by selected cluster key...")

sc.pl.matrixplot(
    rna,
    marker_dict_available,
    groupby=cluster_key,
    use_raw=use_raw,
    standard_scale=dotplot_standard_scale,
    dendrogram=False,
    show=False,
)  # create matrixplot without immediate display so we can save first

save_and_optionally_show_figure(
    output_path=marker_matrixplot_out,
    dpi=dpi,
    display_figures=display_figures,
)  # save and optionally show matrixplot

plt.close()  # release matplotlib figure

print("\nSaved marker matrixplot to:")
print(marker_matrixplot_out)


# ============================================================
# 15. Marker UMAP batches
# ============================================================

print("\nCreating marker UMAP batches...")

marker_umap_files = save_marker_umap_batches(
    adata=rna,
    available_markers=available_markers,
    marker_umaps_dir=marker_umaps_dir,
    batch_size=umap_marker_batch_size,
    use_raw=use_raw,
    dpi=dpi,
    cluster_key=cluster_key_file_token,
    display_figures=display_figures,
)  # save and optionally show marker expression UMAP panels

print("\nSaved marker UMAP files:")
for file_path in marker_umap_files:
    print(file_path)  # terminal audit


# ============================================================
# 16. Marker expression evidence tables
# ============================================================

print("\nComputing mean marker expression by cluster...")

mean_expr = compute_mean_marker_expression(
    adata=rna,
    marker_genes=available_markers,
    cluster_key=cluster_key,
    use_raw=use_raw,
)  # mean expression table

write_csv(
    mean_expr,
    mean_expr_out,
    index=False,
)  # save mean marker expression

print("\nSaved mean marker expression table to:")
print(mean_expr_out)

print("\nComputing percent marker expressing by cluster...")

pct_expr = compute_pct_marker_expressing(
    adata=rna,
    marker_genes=available_markers,
    cluster_key=cluster_key,
    use_raw=use_raw,
)  # percent expressing table

write_csv(
    pct_expr,
    pct_expr_out,
    index=False,
)  # save percent expressing table

print("\nSaved percent-expressing marker table to:")
print(pct_expr_out)


# ============================================================
# 17. Marker group scores and annotation template
# ============================================================

print("\nComputing marker-group scores...")

marker_group_scores = compute_marker_group_scores(
    mean_expr=mean_expr,
    marker_dict_available=marker_dict_available,
)  # simple group-level marker signal per cluster

write_csv(
    marker_group_scores,
    marker_group_scores_out,
    index=False,
)  # save marker group score table

print("\nSaved marker-group score table to:")
print(marker_group_scores_out)

print("\nBuilding manual annotation template...")

annotation_template = build_annotation_template(
    cluster_sizes=cluster_sizes,
    marker_group_scores=marker_group_scores,
)  # worksheet for manual/informed annotation

write_csv(
    annotation_template,
    annotation_template_out,
    index=False,
)  # save manual annotation template

print("\nSaved annotation template table to:")
print(annotation_template_out)

print("\nAnnotation template preview:")
print(annotation_template.head())  # terminal preview


# ============================================================
# 18. Manifest
# ============================================================

manifest = {
    "step": "03A",
    "script": "run_03A_rna_informed_annotation.py",
    "description": (
        "RNA informed marker-based annotation evidence. "
        "This step creates marker evidence outputs only. "
        "It does not write final labels or final annotated h5ad objects."
    ),
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "project": paths.as_dict(),
    "configs": {
        "marker_config": str(marker_config_path),
        "preprocessing_config": str(preprocessing_config_path),
    },
    "clustering_source": {
        "source_config": str(preprocessing_config_path),
        "source_key": "dimensionality_reduction.leiden_key",
        "cluster_key": cluster_key,
        "n_pcs": selected_n_pcs,
        "n_neighbors": selected_n_neighbors,
        "leiden_resolution": selected_leiden_resolution,
    },
    "inputs": {
        "rna_clustered_h5ad": str(rna_clustered_path),
    },
    "outputs": {
        "marker_availability_table": str(marker_availability_out),
        "cluster_sizes_table": str(cluster_sizes_out),
        "marker_dotplot": str(marker_dotplot_out),
        "marker_matrixplot": str(marker_matrixplot_out),
        "marker_umaps_dir": str(marker_umaps_dir),
        "marker_umap_files": marker_umap_files,
        "marker_mean_expression_table": str(mean_expr_out),
        "marker_pct_expressing_table": str(pct_expr_out),
        "marker_group_scores_table": str(marker_group_scores_out),
        "annotation_template_table": str(annotation_template_out),
        "manifest": str(manifest_out),
    },
    "validation": {
        "cluster_key_from_preprocessing_yaml": cluster_key,
        "cluster_key_found_in_rna_obs": bool(cluster_key in rna.obs.columns),
        "umap_found_in_rna_obsm": bool("X_umap" in rna.obsm),
        "all_paths_project_local": True,
        "requested_use_raw": bool(use_raw_requested),
        "actual_use_raw": bool(use_raw),
        "display_figures": bool(display_figures),
        "has_raw": bool(rna.raw is not None),
        "preprocessing_rerun": False,
        "celltypist_run": False,
        "final_annotation_labels_written": False,
        "final_annotated_h5ad_written": False,
    },
    "summary": {
        "n_cells": int(rna.n_obs),
        "n_genes_current_matrix": int(rna.n_vars),
        "n_clusters": int(cluster_sizes.shape[0]),
        "clusters": cluster_sizes["cluster"].astype(str).tolist(),
        "n_marker_groups_requested": int(len(marker_dict)),
        "n_marker_groups_with_available_markers": int(len(marker_dict_available)),
        "marker_groups_requested": list(marker_dict.keys()),
        "marker_groups_with_available_markers": list(marker_dict_available.keys()),
        "n_marker_genes_requested_unique": int(len(marker_genes)),
        "n_marker_genes_available_unique": int(len(available_markers)),
        "n_marker_genes_missing_unique": int(len(missing_markers)),
        "available_markers": available_markers,
        "missing_markers": missing_markers,
    },
    "interpretation_note": (
        "Use combinations of markers rather than single genes. "
        "The annotation template is an evidence worksheet, not final truth. "
        "Mixed or broad clusters should be flagged for review/subclustering. "
        "CellTypist should be used later as complementary evidence."
    ),
}  # complete reproducibility manifest

write_json(
    manifest,
    manifest_out,
)  # save Step 03A manifest

print("\nSaved Step 03A manifest to:")
print(manifest_out)


# ============================================================
# 19. Final summary
# ============================================================

print("\nDone.")
print("Step 03A RNA informed marker-based annotation evidence completed.")
print(f"Selected cluster key used: {cluster_key}")
print("No final labels were added.")
print("No final annotated h5ad was written.")
print("Next: inspect marker evidence and manually review rna_informed_annotation_template.csv.") 



