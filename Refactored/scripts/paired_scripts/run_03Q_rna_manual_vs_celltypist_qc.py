# -*- coding: utf-8 -*-
"""
Step 03Q — Quick RNA manual-marker annotation vs CellTypist QC.

Purpose
-------
- Load the final clustered RNA object from Step 02B.
- Read the selected Leiden key dynamically from preprocessing_rna.yaml.
- Read the manually reviewed labels from the Step 03A annotation template.
- Run CellTypist as independent annotation evidence.
- Compare manual marker-based labels with CellTypist labels.
- Save:
    * RNA object with both annotation columns
    * per-cell annotation table
    * cluster-level comparison table
    * CellTypist label distributions per Leiden cluster
    * selected RNA cluster UMAP
    * manual-label UMAP
    * CellTypist UMAP
    * combined three-panel UMAP
    * manual-vs-CellTypist confusion heatmap

Important
---------
- Manual labels and CellTypist are QC evidence.
- This script does not create final thesis labels automatically.
- Disagreements must be reviewed using marker dotplots and expression tables.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


# ============================================================
# 0. Locate Refactored root robustly
# ============================================================

def find_refactored_root(start_path: Path) -> Path:
    """Find Refactored root containing src/ and projects/."""

    start_path = start_path.resolve()

    for candidate in [start_path, *start_path.parents]:
        if (
            (candidate / "src").is_dir()
            and (candidate / "projects").is_dir()
        ):
            return candidate

    raise FileNotFoundError(
        "Could not locate Refactored root.\n"
        f"Search started from:\n{start_path}"
    )


try:
    script_path = Path(__file__).resolve()
    refactor_dir = find_refactored_root(script_path.parent)
except NameError:
    script_path = Path("run_03Q_rna_manual_vs_celltypist_qc.py")
    refactor_dir = find_refactored_root(Path.cwd())


src_dir = refactor_dir / "src"

if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))


# ============================================================
# 1. Imports from project package
# ============================================================

from pbmcgrn.config import load_yaml
from pbmcgrn.io import write_csv, write_json
from pbmcgrn.paths import ProjectPaths


# ============================================================
# 2. General helper functions
# ============================================================

def assert_path_inside_project(
    path: Path,
    project_root: Path,
    label: str,
) -> None:
    """Fail if a path is outside the active project."""

    resolved_path = Path(path).resolve()
    resolved_root = Path(project_root).resolve()

    try:
        resolved_path.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(
            f"{label} is outside the active project.\n"
            f"Path: {resolved_path}\n"
            f"Project root: {resolved_root}"
        ) from exc


def ensure_parent(path: Path) -> None:
    """Create parent directory for an output file."""

    Path(path).parent.mkdir(
        parents=True,
        exist_ok=True,
    )


def require_section(
    config: dict[str, Any],
    section_name: str,
) -> dict[str, Any]:
    """Return a required YAML dictionary section."""

    if section_name not in config:
        raise KeyError(
            f"Missing YAML section: {section_name}"
        )

    section = config[section_name]

    if not isinstance(section, dict):
        raise TypeError(
            f"YAML section '{section_name}' must be a dictionary."
        )

    return section


def save_show_close(
    path: Path,
    dpi: int = 300,
    show: bool = True,
) -> None:
    """Save current matplotlib figure, optionally show, then close."""

    ensure_parent(path)

    plt.savefig(
        path,
        dpi=dpi,
        bbox_inches="tight",
    )

    if show:
        plt.show()

    plt.close("all")


def normalize_text_label(value: Any) -> str:
    """Convert manual labels to clean strings."""

    if pd.isna(value):
        return ""

    return str(value).strip()


def sort_cluster_table(
    table: pd.DataFrame,
    cluster_column: str,
) -> pd.DataFrame:
    """Sort numeric cluster labels naturally."""

    output = table.copy()

    output["_cluster_numeric"] = pd.to_numeric(
        output[cluster_column],
        errors="coerce",
    )

    output = (
        output
        .sort_values(
            ["_cluster_numeric", cluster_column],
            na_position="last",
        )
        .drop(columns="_cluster_numeric")
        .reset_index(drop=True)
    )

    return output


# ============================================================
# 3. Manual annotation table helpers
# ============================================================

def load_manual_annotation_table(
    annotation_path: Path,
    clusters_in_object: set[str],
) -> pd.DataFrame:
    """
    Load Step 03A annotation template.

    Expected minimum columns:
    - cluster
    - candidate_cell_type

    Optional:
    - confidence
    - notes
    - needs_subclustering
    """

    if not annotation_path.exists():
        raise FileNotFoundError(
            "Step 03A annotation template was not found:\n"
            f"{annotation_path}"
        )

    annotation = pd.read_csv(
    annotation_path,
    sep=None,
    engine="python",
    encoding="utf-8-sig",
    dtype=str,
    )

    annotation.columns = (
    annotation.columns
    .astype(str)
    .str.replace("\ufeff", "", regex=False)
    .str.strip()
    )

    print("\nDetected annotation columns:")
    print(annotation.columns.tolist())
    required_columns = {
        "cluster",
        "candidate_cell_type_broad",
    }

    missing_columns = required_columns.difference(
        annotation.columns
    )

    if missing_columns:
        raise KeyError(
            "Manual annotation table is missing columns:\n"
            f"{sorted(missing_columns)}"
        )

    annotation["cluster"] = (
        annotation["cluster"]
        .astype(str)
        .str.strip()
    )

    annotation["candidate_cell_type_broad"] = (
        annotation["candidate_cell_type_broad"]
        .apply(normalize_text_label)
)

    if annotation["cluster"].duplicated().any():
        duplicated = (
            annotation.loc[
                annotation["cluster"].duplicated(),
                "cluster",
            ]
            .tolist()
        )

        raise ValueError(
            "Manual annotation table has duplicated clusters:\n"
            f"{duplicated}"
        )

    missing_labels = annotation.loc[
        annotation["candidate_cell_type_broad"] == "",
        "cluster",
    ].tolist()

    if missing_labels:
        raise ValueError(
            "Some clusters do not yet have candidate_cell_type labels.\n"
            f"Unlabeled clusters: {missing_labels}\n\n"
            "Fill the candidate_cell_type column in:\n"
            f"{annotation_path}"
        )

    clusters_in_table = set(
        annotation["cluster"].astype(str)
    )

    missing_from_table = sorted(
        clusters_in_object - clusters_in_table
    )

    extra_in_table = sorted(
        clusters_in_table - clusters_in_object
    )

    if missing_from_table:
        raise ValueError(
            "Clusters present in RNA but missing from manual table:\n"
            f"{missing_from_table}"
        )

    if extra_in_table:
        raise ValueError(
            "Clusters present in manual table but absent from RNA:\n"
            f"{extra_in_table}"
        )

    optional_defaults = {
        "confidence": "",
        "notes": "",
        "needs_subclustering": "",
    }

    for column, default_value in optional_defaults.items():
        if column not in annotation.columns:
            annotation[column] = default_value

    return sort_cluster_table(
        annotation,
        cluster_column="cluster",
    )


# ============================================================
# 4. CellTypist helpers
# ============================================================

def prepare_celltypist_input(
    rna: sc.AnnData,
) -> sc.AnnData:
    """
    Prepare log-normalized full-gene RNA for CellTypist.

    The Step 02B object stores the complete log-normalized
    gene matrix in rna.raw.
    """

    if rna.raw is not None:
        celltypist_input = rna.raw.to_adata()
    else:
        celltypist_input = rna.copy()

    # Preserve the original UMAP and cluster metadata.
    celltypist_input.obs = rna.obs.copy()

    if "X_umap" in rna.obsm:
        celltypist_input.obsm["X_umap"] = (
            rna.obsm["X_umap"].copy()
        )

    if not celltypist_input.obs_names.equals(
        rna.obs_names
    ):
        raise RuntimeError(
            "CellTypist input does not retain the same cell order."
        )

    return celltypist_input


def run_celltypist_annotation(
    rna: sc.AnnData,
    model_name: str,
) -> pd.DataFrame:
    """Run CellTypist and return aligned per-cell predictions."""

    try:
        import celltypist
        from celltypist import models
    except ImportError as exc:
        raise ImportError(
            "CellTypist is not installed.\n"
            "Install it in the current environment with:\n"
            "pip install celltypist"
        ) from exc

    print("\nCellTypist version:")
    print(celltypist.__version__)

    print("\nChecking CellTypist model availability...")

    models.download_models(
        model=model_name,
        force_update=False,
    )

    celltypist_input = prepare_celltypist_input(
        rna
    )

    print("\nCellTypist input:")
    print(celltypist_input)

    print("\nRunning CellTypist model:")
    print(model_name)

    predictions = celltypist.annotate(
        celltypist_input,
        model=model_name,
        majority_voting=True,
    )

    prediction_adata = predictions.to_adata()

    prediction_obs = prediction_adata.obs.copy()

    if not prediction_obs.index.equals(
        rna.obs_names
    ):
        prediction_obs = prediction_obs.reindex(
            rna.obs_names
        )

    if prediction_obs.index.hasnans:
        raise RuntimeError(
            "CellTypist predictions could not be aligned "
            "to the RNA cell barcodes."
        )

    if "predicted_labels" not in prediction_obs.columns:
        raise KeyError(
            "CellTypist output does not contain predicted_labels.\n"
            f"Available columns: {prediction_obs.columns.tolist()}"
        )

    output = pd.DataFrame(
        index=rna.obs_names
    )

    output["celltypist_predicted_label"] = (
        prediction_obs["predicted_labels"]
        .astype(str)
        .values
    )

    if "majority_voting" in prediction_obs.columns:
        output["celltypist_majority_voting"] = (
            prediction_obs["majority_voting"]
            .astype(str)
            .values
        )
    else:
        output["celltypist_majority_voting"] = (
            output["celltypist_predicted_label"]
        )

    confidence_candidates = [
        "conf_score",
        "probability",
        "score",
    ]

    confidence_column = next(
        (
            column
            for column in confidence_candidates
            if column in prediction_obs.columns
        ),
        None,
    )

    if confidence_column is not None:
        output["celltypist_confidence"] = (
            prediction_obs[confidence_column]
            .values
        )
    else:
        output["celltypist_confidence"] = pd.NA

    return output


# ============================================================
# 5. Comparison helpers
# ============================================================

def build_celltypist_cluster_summary(
    rna: sc.AnnData,
    cluster_key: str,
) -> pd.DataFrame:
    """Summarize CellTypist label composition per Leiden cluster."""

    summary = (
        rna.obs
        .groupby(
            [
                cluster_key,
                "celltypist_majority_voting",
            ],
            observed=False,
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
    )

    totals = (
        rna.obs
        .groupby(
            cluster_key,
            observed=False,
            dropna=False,
        )
        .size()
        .rename("cluster_n_cells")
        .reset_index()
    )

    summary = summary.merge(
        totals,
        on=cluster_key,
        how="left",
    )

    summary["fraction_in_cluster"] = (
        summary["n_cells"]
        / summary["cluster_n_cells"]
    )

    summary[cluster_key] = (
        summary[cluster_key]
        .astype(str)
    )

    summary = (
        summary
        .sort_values(
            [
                cluster_key,
                "n_cells",
            ],
            ascending=[True, False],
        )
        .reset_index(drop=True)
    )

    return sort_cluster_table(
        summary,
        cluster_column=cluster_key,
    )


def get_celltypist_majority_per_cluster(
    cluster_summary: pd.DataFrame,
    cluster_key: str,
) -> pd.DataFrame:
    """Select the most frequent CellTypist label in every cluster."""

    majority = (
        cluster_summary
        .sort_values(
            [
                cluster_key,
                "n_cells",
            ],
            ascending=[True, False],
        )
        .groupby(
            cluster_key,
            as_index=False,
            sort=False,
        )
        .head(1)
        .copy()
    )

    majority = majority.rename(
        columns={
            "celltypist_majority_voting":
                "celltypist_cluster_majority_label",
            "n_cells":
                "celltypist_cluster_majority_n_cells",
            "fraction_in_cluster":
                "celltypist_cluster_majority_fraction",
        }
    )

    return sort_cluster_table(
        majority,
        cluster_column=cluster_key,
    )


def classify_annotation_agreement(
    manual_label: str,
    celltypist_label: str,
) -> str:
    """
    Conservative string-level agreement classification.

    Exact biological harmonization should be reviewed manually.
    """

    manual = str(manual_label).strip().lower()
    automatic = str(celltypist_label).strip().lower()

    if manual == automatic:
        return "exact"

    manual_tokens = set(
        manual
        .replace("-", " ")
        .replace("_", " ")
        .split()
    )

    automatic_tokens = set(
        automatic
        .replace("-", " ")
        .replace("_", " ")
        .split()
    )

    shared_tokens = manual_tokens.intersection(
        automatic_tokens
    )

    biologically_informative_tokens = {
        "cd4",
        "cd8",
        "t",
        "b",
        "nk",
        "monocyte",
        "monocytes",
        "dendritic",
        "platelet",
        "plasma",
        "naive",
        "memory",
        "cytotoxic",
    }

    if shared_tokens.intersection(
        biologically_informative_tokens
    ):
        return "partial"

    return "different"


def build_cluster_comparison(
    manual_annotation: pd.DataFrame,
    celltypist_majority: pd.DataFrame,
    cluster_key: str,
) -> pd.DataFrame:
    """Merge manual and CellTypist cluster labels."""

    manual = manual_annotation.copy()

    manual = manual.rename(
    columns={
        "cluster": cluster_key,
        "candidate_cell_type_broad":
            "manual_marker_label",
        "annotation_confidence":
            "manual_annotation_confidence",
        "review_notes":
            "manual_annotation_notes",
        }
    )

    manual[cluster_key] = (
        manual[cluster_key].astype(str)
    )

    celltypist_majority[cluster_key] = (
        celltypist_majority[cluster_key]
        .astype(str)
    )

    comparison = manual.merge(
        celltypist_majority[
            [
                cluster_key,
                "celltypist_cluster_majority_label",
                "celltypist_cluster_majority_n_cells",
                "cluster_n_cells",
                "celltypist_cluster_majority_fraction",
            ]
        ],
        on=cluster_key,
        how="left",
    )

    comparison["annotation_agreement"] = [
        classify_annotation_agreement(
            manual_label,
            automatic_label,
        )
        for manual_label, automatic_label in zip(
            comparison["manual_marker_label"],
            comparison[
                "celltypist_cluster_majority_label"
            ],
        )
    ]

    comparison["manual_review_needed"] = (
        (
            comparison["annotation_agreement"]
            != "exact"
        )
        |
        (
            comparison[
                "celltypist_cluster_majority_fraction"
            ]
            < 0.75
        )
    )

    return sort_cluster_table(
        comparison,
        cluster_column=cluster_key,
    )


# ============================================================
# 6. Main workflow
# ============================================================

def main() -> None:
    # --------------------------------------------------------
    # 6.1 Active project
    # --------------------------------------------------------

    paths = ProjectPaths(
        refactor_root=refactor_dir
    )

    paths.ensure_output_dirs()

    expected_project_id = (
        "pbmc_10k_multiome_v1"
    )

    print("\n" + "=" * 80)
    print("STEP 03Q — MANUAL MARKER LABELS VS CELLTYPIST QC")
    print("=" * 80)

    print("\nRefactored root:")
    print(refactor_dir)

    print("\nActive project:")
    print(paths.project_id)

    print("\nProject root:")
    print(paths.project_root)

    if paths.project_id != expected_project_id:
        raise RuntimeError(
            "Wrong active project.\n"
            f"Expected: {expected_project_id}\n"
            f"Observed: {paths.project_id}"
        )

    # --------------------------------------------------------
    # 6.2 Read final clustering source
    # --------------------------------------------------------

    preprocessing_path = (
        paths.configs
        / "preprocessing_rna.yaml"
    )

    if not preprocessing_path.exists():
        raise FileNotFoundError(
            "preprocessing_rna.yaml not found:\n"
            f"{preprocessing_path}"
        )

    preprocessing_config = load_yaml(
        preprocessing_path
    )

    dim_cfg = require_section(
        preprocessing_config,
        "dimensionality_reduction",
    )

    cluster_key = str(
        dim_cfg["leiden_key"]
    )

    print("\nSelected clustering:")
    print(f"cluster_key: {cluster_key}")
    print(f"n_pcs: {dim_cfg['n_pcs']}")
    print(f"n_neighbors: {dim_cfg['n_neighbors']}")
    print(
        "Leiden resolution:",
        dim_cfg["leiden_resolution"],
    )

    # --------------------------------------------------------
    # 6.3 Input paths
    # --------------------------------------------------------

    rna_clustered_path = (
        paths.processed_rna
        / "pbmc_multiome_rna_clustered.h5ad"
    )

    manual_annotation_path = (
        paths.project_root
        / "reports"
        / "tables"
        / "rna"
        / "step03A"
        / "rna_informed_annotation_template.csv"
    )

    for label, path in {
        "rna_clustered": rna_clustered_path,
        "manual_annotation": manual_annotation_path,
        "preprocessing_config": preprocessing_path,
    }.items():
        assert_path_inside_project(
            path,
            paths.project_root,
            label,
        )

        if not path.exists():
            raise FileNotFoundError(
                f"Required input not found ({label}):\n{path}"
            )

    # --------------------------------------------------------
    # 6.4 Output paths
    # --------------------------------------------------------

    output_table_dir = (
        paths.project_root
        / "reports"
        / "tables"
        / "rna"
        / "step03Q"
    )

    output_figure_dir = (
        paths.project_root
        / "reports"
        / "figures"
        / "rna"
        / "step03Q"
    )

    output_table_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    annotated_h5ad_out = (
        paths.processed_rna
        / "pbmc_multiome_rna_manual_celltypist_qc.h5ad"
    )

    per_cell_out = (
        output_table_dir
        / "rna_manual_celltypist_per_cell.csv"
    )

    celltypist_cluster_summary_out = (
        output_table_dir
        / "rna_celltypist_label_distribution_by_cluster.csv"
    )

    celltypist_majority_out = (
        output_table_dir
        / "rna_celltypist_majority_label_by_cluster.csv"
    )

    cluster_comparison_out = (
        output_table_dir
        / "rna_manual_vs_celltypist_cluster_comparison.csv"
    )

    selected_clusters_umap_out = (
        output_figure_dir
        / f"rna_selected_clusters_{cluster_key}.png"
    )

    manual_umap_out = (
        output_figure_dir
        / f"rna_manual_marker_labels_{cluster_key}.png"
    )

    celltypist_umap_out = (
        output_figure_dir
        / f"rna_celltypist_labels_{cluster_key}.png"
    )

    combined_umap_out = (
        output_figure_dir
        / f"rna_selected_clusters_manual_celltypist_{cluster_key}.png"
    )

    confusion_heatmap_out = (
        output_figure_dir
        / f"rna_manual_vs_celltypist_heatmap_{cluster_key}.png"
    )

    manifest_out = (
        paths.project_root
        / "reports"
        / "manifests"
        / "step03Q_rna_manual_vs_celltypist_qc_manifest.json"
    )

    output_paths = {
        "annotated_h5ad": annotated_h5ad_out,
        "per_cell_table": per_cell_out,
        "celltypist_cluster_summary":
            celltypist_cluster_summary_out,
        "celltypist_majority":
            celltypist_majority_out,
        "cluster_comparison":
            cluster_comparison_out,
        "selected_clusters_umap":
            selected_clusters_umap_out,
        "manual_labels_umap":
            manual_umap_out,
        "celltypist_umap":
            celltypist_umap_out,
        "combined_umap":
            combined_umap_out,
        "confusion_heatmap":
            confusion_heatmap_out,
        "manifest":
            manifest_out,
    }

    for label, path in output_paths.items():
        assert_path_inside_project(
            path,
            paths.project_root,
            label,
        )

        ensure_parent(path)

    # --------------------------------------------------------
    # 6.5 Load clustered RNA
    # --------------------------------------------------------

    print("\nLoading clustered RNA:")
    print(rna_clustered_path)

    rna = sc.read_h5ad(
        rna_clustered_path
    )

    print("\nRNA object:")
    print(rna)

    if cluster_key not in rna.obs.columns:
        raise KeyError(
            f"Cluster key not found in rna.obs: {cluster_key}"
        )

    if "X_umap" not in rna.obsm:
        raise KeyError(
            "X_umap is missing from the clustered RNA object."
        )

    rna.obs[cluster_key] = (
        rna.obs[cluster_key]
        .astype(str)
        .astype("category")
    )

    clusters_in_object = set(
        rna.obs[cluster_key]
        .astype(str)
        .unique()
    )

    print("\nSelected RNA cluster counts:")
    print(
        rna.obs[cluster_key]
        .value_counts()
        .sort_index()
    )

    # --------------------------------------------------------
    # 6.6 Load manual marker labels
    # --------------------------------------------------------

    manual_annotation = (
        load_manual_annotation_table(
            annotation_path=manual_annotation_path,
            clusters_in_object=clusters_in_object,
        )
    )

    manual_mapping = (
        manual_annotation
        .set_index("cluster")
        ["candidate_cell_type_broad"]
        .to_dict()
        )
    confidence_mapping = (
        manual_annotation
        .set_index("cluster")
        ["confidence"]
        .to_dict()
    )

    notes_mapping = (
        manual_annotation
        .set_index("cluster")
        ["notes"]
        .to_dict()
    )

    rna.obs["manual_marker_label"] = (
        rna.obs[cluster_key]
        .astype(str)
        .map(manual_mapping)
        .astype("category")
    )

    rna.obs["manual_annotation_confidence"] = (
        rna.obs[cluster_key]
        .astype(str)
        .map(confidence_mapping)
    )

    rna.obs["manual_annotation_notes"] = (
        rna.obs[cluster_key]
        .astype(str)
        .map(notes_mapping)
    )

    if rna.obs["manual_marker_label"].isna().any():
        raise RuntimeError(
            "Some RNA cells did not receive a manual marker label."
        )

    print("\nManual marker labels:")
    print(
        rna.obs["manual_marker_label"]
        .value_counts()
    )

    # --------------------------------------------------------
    # 6.7 Run CellTypist
    # --------------------------------------------------------

    celltypist_model = (
        "Immune_All_Low.pkl"
    )

    celltypist_results = (
        run_celltypist_annotation(
            rna=rna,
            model_name=celltypist_model,
        )
    )

    for column in celltypist_results.columns:
        rna.obs[column] = (
            celltypist_results[column]
            .reindex(rna.obs_names)
            .values
        )

    rna.obs["celltypist_majority_voting"] = (
        rna.obs["celltypist_majority_voting"]
        .astype(str)
        .astype("category")
    )

    print("\nCellTypist majority-voting labels:")
    print(
        rna.obs["celltypist_majority_voting"]
        .value_counts()
    )

    # --------------------------------------------------------
    # 6.8 Build comparison tables
    # --------------------------------------------------------

    celltypist_cluster_summary = (
        build_celltypist_cluster_summary(
            rna=rna,
            cluster_key=cluster_key,
        )
    )

    celltypist_majority = (
        get_celltypist_majority_per_cluster(
            cluster_summary=celltypist_cluster_summary,
            cluster_key=cluster_key,
        )
    )

    cluster_comparison = (
        build_cluster_comparison(
            manual_annotation=manual_annotation,
            celltypist_majority=celltypist_majority,
            cluster_key=cluster_key,
        )
    )

    write_csv(
        celltypist_cluster_summary,
        celltypist_cluster_summary_out,
        index=False,
    )

    write_csv(
        celltypist_majority,
        celltypist_majority_out,
        index=False,
    )

    write_csv(
        cluster_comparison,
        cluster_comparison_out,
        index=False,
    )

    print("\nManual vs CellTypist comparison:")
    print(
        cluster_comparison[
            [
                cluster_key,
                "manual_marker_label",
                "celltypist_cluster_majority_label",
                "celltypist_cluster_majority_fraction",
                "annotation_agreement",
                "manual_review_needed",
            ]
        ].to_string(index=False)
    )

    # --------------------------------------------------------
    # 6.9 Save per-cell table
    # --------------------------------------------------------

    per_cell = rna.obs[
        [
            cluster_key,
            "manual_marker_label",
            "manual_annotation_confidence",
            "celltypist_predicted_label",
            "celltypist_majority_voting",
            "celltypist_confidence",
        ]
    ].copy()

    per_cell.insert(
        0,
        "cell_barcode",
        per_cell.index.astype(str),
    )

    write_csv(
        per_cell,
        per_cell_out,
        index=False,
    )

    # --------------------------------------------------------
    # 6.10 Save RNA with both label systems
    # --------------------------------------------------------

    rna.write_h5ad(
        annotated_h5ad_out
    )

    print("\nSaved RNA object with manual and CellTypist labels:")
    print(annotated_h5ad_out)

    # --------------------------------------------------------
    # 6.11 UMAP: selected Leiden clusters
    # --------------------------------------------------------

    sc.pl.umap(
        rna,
        color=cluster_key,
        legend_loc="on data",
        legend_fontsize=7,
        legend_fontweight="normal",
        frameon=False,
        title="Selected RNA Leiden clusters",
        show=False,
    )

    save_show_close(
        selected_clusters_umap_out,
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.12 UMAP: manual marker labels
    # --------------------------------------------------------

    sc.pl.umap(
        rna,
        color="manual_marker_label",
        legend_loc="right margin",
        frameon=False,
        title="Manual marker-based annotation",
        show=False,
    )

    save_show_close(
        manual_umap_out,
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.13 UMAP: CellTypist labels
    # --------------------------------------------------------

    sc.pl.umap(
        rna,
        color="celltypist_majority_voting",
        legend_loc="right margin",
        frameon=False,
        title="CellTypist majority-voting annotation",
        show=False,
    )

    save_show_close(
        celltypist_umap_out,
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.14 Combined UMAP panel requested by user
    # --------------------------------------------------------

    sc.pl.umap(
        rna,
        color=[
            cluster_key,
            "manual_marker_label",
            "celltypist_majority_voting",
        ],
        legend_loc="right margin",
        frameon=False,
        title=[
            "Selected RNA clusters",
            "Manual marker labeling",
            "CellTypist labeling",
        ],
        ncols=3,
        show=False,
    )

    save_show_close(
        combined_umap_out,
        dpi=300,
        show=True,
    )

    print("\nSaved combined selected-cluster/manual/CellTypist UMAP:")
    print(combined_umap_out)

    # --------------------------------------------------------
    # 6.15 Manual-vs-CellTypist confusion heatmap
    # --------------------------------------------------------

    confusion = pd.crosstab(
        rna.obs["manual_marker_label"].astype(str),
        rna.obs[
            "celltypist_majority_voting"
        ].astype(str),
        normalize="index",
    )

    plt.figure(
        figsize=(
            max(10, confusion.shape[1] * 0.55),
            max(6, confusion.shape[0] * 0.55),
        )
    )

    sns.heatmap(
        confusion,
        cmap="Blues",
        vmin=0,
        vmax=1,
        annot=True,
        fmt=".2f",
        linewidths=0.3,
        linecolor="white",
        cbar_kws={
            "label": (
                "Fraction within manual marker label"
            )
        },
    )

    plt.title(
        "Manual marker annotation vs CellTypist"
    )

    plt.xlabel(
        "CellTypist majority-voting label"
    )

    plt.ylabel(
        "Manual marker-based label"
    )

    plt.xticks(
        rotation=45,
        ha="right",
    )

    plt.yticks(
        rotation=0,
    )

    save_show_close(
        confusion_heatmap_out,
        dpi=300,
        show=True,
    )

    # --------------------------------------------------------
    # 6.16 Manifest
    # --------------------------------------------------------

    manifest = {
        "step": "03Q",
        "description": (
            "Quick RNA annotation QC comparing manual "
            "marker-based labels with CellTypist."
        ),
        "timestamp": datetime.now().isoformat(
            timespec="seconds"
        ),
        "project": paths.as_dict(),
        "inputs": {
            "rna_clustered_h5ad": str(
                rna_clustered_path
            ),
            "manual_annotation_table": str(
                manual_annotation_path
            ),
            "preprocessing_config": str(
                preprocessing_path
            ),
        },
        "clustering": {
            "cluster_key": cluster_key,
            "n_pcs": int(
                dim_cfg["n_pcs"]
            ),
            "n_neighbors": int(
                dim_cfg["n_neighbors"]
            ),
            "leiden_resolution": float(
                dim_cfg["leiden_resolution"]
            ),
        },
        "celltypist": {
            "model": celltypist_model,
            "majority_voting": True,
            "role": (
                "annotation confirmation evidence, "
                "not final biological truth"
            ),
        },
        "summary": {
            "n_cells": int(rna.n_obs),
            "n_clusters": int(
                rna.obs[cluster_key].nunique()
            ),
            "n_manual_labels": int(
                rna.obs[
                    "manual_marker_label"
                ].nunique()
            ),
            "n_celltypist_labels": int(
                rna.obs[
                    "celltypist_majority_voting"
                ].nunique()
            ),
            "clusters_requiring_review": int(
                cluster_comparison[
                    "manual_review_needed"
                ].sum()
            ),
        },
        "outputs": {
            label: str(path)
            for label, path in output_paths.items()
        },
    }

    write_json(
        manifest,
        manifest_out,
    )

    # --------------------------------------------------------
    # 6.17 Final report
    # --------------------------------------------------------

    print("\n" + "=" * 80)
    print("STEP 03Q COMPLETED")
    print("=" * 80)

    print("\nCombined UMAP:")
    print(combined_umap_out)

    print("\nCluster comparison table:")
    print(cluster_comparison_out)

    print("\nConfusion heatmap:")
    print(confusion_heatmap_out)

    print("\nAnnotated RNA object:")
    print(annotated_h5ad_out)

    print(
        "\nInterpret disagreements together with the Step 03A "
        "dotplot, matrixplot, marker UMAPs, mean-expression "
        "and percent-expressing tables."
    )


if __name__ == "__main__":
    main()