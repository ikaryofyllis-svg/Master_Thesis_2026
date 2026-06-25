"""Βήμα 03B — Sanity check annotation strategies, Leiden-reference agreement και group decision table"""

from pathlib import Path
import sys
import json

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# ============================================================
# 0. Refactored project root και Python import path
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


paths = ProjectPaths(refactor_root=refactor_dir)
paths.ensure_output_dirs()

print("Refactored project root:")
print(paths.refactor_root)


# ============================================================
# 2. Parameters για annotation sanity check
# ============================================================

cluster_key = "leiden_0_6"

high_consensus_threshold = 0.80
medium_consensus_threshold = 0.60
split_secondary_fraction_threshold = 0.25
small_cluster_n_cells = 100
top_de_n_genes = 10

print("\nDecision thresholds:")
print("High consensus threshold:", high_consensus_threshold)
print("Medium consensus threshold:", medium_consensus_threshold)
print("Split secondary fraction threshold:", split_secondary_fraction_threshold)
print("Small cluster threshold:", small_cluster_n_cells)


# ============================================================
# 3. Helper function για resolving paths
# ============================================================

def resolve_from_refactor(path_string):
    """
    Resolve relative paths from Refactored/.
    """

    path = Path(path_string)

    if path.is_absolute():
        return path

    return (paths.refactor_root / path).resolve()


# ============================================================
# 4. Find reference-annotated RNA object
# ============================================================
# Υποστηρίζει και το νέο generic output και το προηγούμενο CellTypist-specific output.
#
# Preferred:
#   pbmc_rna_reference_annotated.h5ad
#
# Fallback:
#   pbmc_rna_celltypist_annotated.h5ad

generic_reference_h5ad = paths.processed_rna / "pbmc_rna_reference_annotated.h5ad"
celltypist_specific_h5ad = paths.processed_rna / "pbmc_rna_celltypist_annotated.h5ad"

if generic_reference_h5ad.exists():
    reference_h5ad = generic_reference_h5ad
    reference_source_name = "generic_reference_backend"

elif celltypist_specific_h5ad.exists():
    reference_h5ad = celltypist_specific_h5ad
    reference_source_name = "celltypist_specific_backend"

else:
    raise FileNotFoundError(
        "Δεν βρέθηκε reference-annotated RNA object.\n"
        "Πρέπει πρώτα να τρέξει ένα από:\n"
        f"{generic_reference_h5ad}\n"
        f"{celltypist_specific_h5ad}"
    )

print("\nLoading reference-annotated object:")
print(reference_h5ad)
print("Reference source:", reference_source_name)

adata = sc.read_h5ad(reference_h5ad)

print("\nReference-annotated object:")
print(adata)

print("\nobs columns:")
print(adata.obs.columns.tolist())


# ============================================================
# 5. Standardize reference label column
# ============================================================
# Για generic backend περιμένουμε:
#   reference_label
#
# Για CellTypist-specific παλιό script:
#   majority_voting ή predicted_labels
#
# Εδώ φτιάχνουμε πάντα:
#   adata.obs["reference_label"]

if cluster_key not in adata.obs.columns:
    raise ValueError(f"Missing Leiden cluster key: {cluster_key}")

if "reference_label" in adata.obs.columns:
    reference_label_col = "reference_label"

elif "majority_voting" in adata.obs.columns:
    reference_label_col = "majority_voting"
    adata.obs["reference_label"] = adata.obs["majority_voting"].astype(str)

elif "predicted_labels" in adata.obs.columns:
    reference_label_col = "predicted_labels"
    adata.obs["reference_label"] = adata.obs["predicted_labels"].astype(str)

else:
    raise ValueError(
        "No usable reference label column found.\n"
        "Expected one of: reference_label, majority_voting, predicted_labels."
    )

adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
adata.obs["reference_label"] = adata.obs["reference_label"].astype(str)

print("\nUsing reference label column:")
print(reference_label_col)

print("\nNumber of Leiden clusters:")
print(adata.obs[cluster_key].nunique())

print("\nNumber of reference labels:")
print(adata.obs["reference_label"].nunique())


# ============================================================
# 6. Strategy-level group counts
# ============================================================

leiden_counts = (
    adata.obs[cluster_key]
    .value_counts()
    .rename_axis("group")
    .reset_index(name="n_cells")
)

leiden_counts.insert(0, "strategy", "leiden_clusters")

reference_counts = (
    adata.obs["reference_label"]
    .value_counts()
    .rename_axis("group")
    .reset_index(name="n_cells")
)

reference_counts.insert(0, "strategy", "reference_labels")

strategy_counts = pd.concat(
    [leiden_counts, reference_counts],
    axis=0,
    ignore_index=True,
)

strategy_counts_out = (
    paths.tables / "rna_annotation_strategy_group_counts.csv"
)

write_csv(strategy_counts, strategy_counts_out, index=False)

print("\nSaved strategy group counts to:")
print(strategy_counts_out)


# ============================================================
# 7. Leiden × reference label counts/fractions
# ============================================================

count_table = pd.crosstab(
    adata.obs[cluster_key],
    adata.obs["reference_label"],
)

count_table.index.name = "leiden_cluster"

fraction_by_leiden = count_table.div(
    count_table.sum(axis=1),
    axis=0,
)

fraction_by_reference = count_table.div(
    count_table.sum(axis=0),
    axis=1,
)

count_table_out = (
    paths.tables / "rna_sanity_leiden_reference_counts.csv"
)

fraction_by_leiden_out = (
    paths.tables / "rna_sanity_leiden_reference_fractions_by_leiden.csv"
)

fraction_by_reference_out = (
    paths.tables / "rna_sanity_reference_fractions_by_reference_label.csv"
)

count_table.to_csv(count_table_out)
fraction_by_leiden.to_csv(fraction_by_leiden_out)
fraction_by_reference.to_csv(fraction_by_reference_out)

print("\nSaved Leiden × reference count/fraction tables.")


# ============================================================
# 8. Per-Leiden decision table
# ============================================================

decision_rows = []

for cluster, row in count_table.iterrows():
    n_cells = int(row.sum())

    sorted_counts = row.sort_values(ascending=False)
    sorted_fractions = sorted_counts / n_cells

    labels = sorted_counts.index.tolist()
    counts = sorted_counts.values.tolist()
    fractions = sorted_fractions.values.tolist()

    top_label = labels[0]
    top_count = int(counts[0])
    top_fraction = float(fractions[0])

    second_label = labels[1] if len(labels) > 1 else ""
    second_count = int(counts[1]) if len(labels) > 1 else 0
    second_fraction = float(fractions[1]) if len(labels) > 1 else 0.0

    n_labels_ge_10pct = int((sorted_fractions >= 0.10).sum())
    n_labels_ge_25pct = int((sorted_fractions >= split_secondary_fraction_threshold).sum())

    if top_fraction >= high_consensus_threshold:
        consensus_status = "high_consensus"
    elif top_fraction >= medium_consensus_threshold:
        consensus_status = "medium_consensus"
    else:
        consensus_status = "mixed_or_low_consensus"

    if n_cells < small_cluster_n_cells:
        size_status = "small_cluster_review"
    else:
        size_status = "adequate_size"

    if consensus_status == "high_consensus":
        recommended_action = "keep_leiden_cluster_assign_top_label"
        use_for_downstream_default = True

    elif consensus_status == "medium_consensus" and second_fraction < split_secondary_fraction_threshold:
        recommended_action = "keep_with_manual_review"
        use_for_downstream_default = True

    elif second_fraction >= split_secondary_fraction_threshold:
        recommended_action = "review_possible_split_or_subcluster"
        use_for_downstream_default = False

    else:
        recommended_action = "ambiguous_review_or_exclude"
        use_for_downstream_default = False

    if size_status == "small_cluster_review":
        recommended_action = recommended_action + "__small_cluster_caution"

    decision_rows.append(
        {
            "leiden_cluster": str(cluster),
            "n_cells": n_cells,
            "top_reference_label": str(top_label),
            "top_reference_count": top_count,
            "top_reference_fraction": top_fraction,
            "second_reference_label": str(second_label),
            "second_reference_count": second_count,
            "second_reference_fraction": second_fraction,
            "n_reference_labels_ge_10pct": n_labels_ge_10pct,
            "n_reference_labels_ge_25pct": n_labels_ge_25pct,
            "top_5_reference_labels": "; ".join(map(str, labels[:5])),
            "top_5_reference_fractions": "; ".join([f"{x:.3f}" for x in fractions[:5]]),
            "consensus_status": consensus_status,
            "size_status": size_status,
            "recommended_action": recommended_action,
            "use_for_downstream_default": use_for_downstream_default,
            "manual_decision": "",
            "manual_final_label": "",
            "manual_notes": "",
            "use_for_downstream_final": "",
        }
    )

decision_table = pd.DataFrame(decision_rows)

decision_table = decision_table.sort_values(
    by="leiden_cluster",
    key=lambda x: x.astype(int),
).reset_index(drop=True)


# ============================================================
# 9. Reference label spread table
# ============================================================
# Αυτό δείχνει αν ένα reference label είναι σπασμένο σε πολλά Leiden clusters.
# Αν ναι, πιθανό downstream merge candidate.

spread_rows = []

for reference_label in count_table.columns:
    col_counts = count_table[reference_label]
    total_cells_for_label = int(col_counts.sum())

    if total_cells_for_label == 0:
        continue

    fractions_across_leiden = col_counts / total_cells_for_label
    nonzero_clusters = col_counts[col_counts > 0]

    top_leiden_cluster = fractions_across_leiden.sort_values(ascending=False).index[0]
    top_leiden_fraction = float(fractions_across_leiden.max())

    leiden_clusters_ge_10pct = (
        fractions_across_leiden[fractions_across_leiden >= 0.10]
        .index
        .astype(str)
        .tolist()
    )

    n_leiden_clusters_ge_10pct = len(leiden_clusters_ge_10pct)

    if n_leiden_clusters_ge_10pct > 1:
        merge_signal = "reference_label_spans_multiple_leiden_clusters"
    else:
        merge_signal = "mostly_one_leiden_cluster"

    spread_rows.append(
        {
            "reference_label": str(reference_label),
            "n_cells_reference_label": total_cells_for_label,
            "n_leiden_clusters_with_any_cells": int(nonzero_clusters.shape[0]),
            "n_leiden_clusters_ge_10pct_of_reference_label": n_leiden_clusters_ge_10pct,
            "leiden_clusters_ge_10pct_of_reference_label": "; ".join(leiden_clusters_ge_10pct),
            "top_leiden_cluster_for_reference_label": str(top_leiden_cluster),
            "top_leiden_fraction_of_reference_label": top_leiden_fraction,
            "merge_signal": merge_signal,
        }
    )

reference_spread = pd.DataFrame(spread_rows)

reference_spread = reference_spread.sort_values(
    by="n_cells_reference_label",
    ascending=False,
).reset_index(drop=True)


# ============================================================
# 10. Compute top DE genes per Leiden cluster
# ============================================================
# Δεν χρησιμοποιούμε χειροκίνητη marker list.
# Τα genes προκύπτουν από το ίδιο το dataset.
#
# Αν το adata είναι full-gene log-normalized, είμαστε καλά.
# Αν είναι CellTypist output, συνήθως έχει log-normalized input από Step 03A.

print("\nComputing top DE genes per Leiden cluster...")

sc.tl.rank_genes_groups(
    adata,
    groupby=cluster_key,
    method="wilcoxon",
    use_raw=False,
)

de_rows = []

cluster_categories = sorted(
    adata.obs[cluster_key].astype(str).unique(),
    key=lambda x: int(x) if x.isdigit() else x,
)

for cluster in cluster_categories:
    de_df = sc.get.rank_genes_groups_df(
        adata,
        group=cluster,
    )

    de_df = de_df.sort_values("scores", ascending=False).head(top_de_n_genes)

    for rank, row in enumerate(de_df.itertuples(index=False), start=1):
        de_rows.append(
            {
                "leiden_cluster": str(cluster),
                "rank": rank,
                "gene": row.names,
                "score": row.scores,
                "logfoldchange": row.logfoldchanges,
                "pvals_adj": row.pvals_adj,
            }
        )

top_de = pd.DataFrame(de_rows)

top_de_compact = (
    top_de
    .groupby("leiden_cluster")["gene"]
    .apply(lambda genes: "; ".join(genes.astype(str).tolist()))
    .reset_index(name=f"top_{top_de_n_genes}_de_genes")
)

decision_table = decision_table.merge(
    top_de_compact,
    on="leiden_cluster",
    how="left",
)


# ============================================================
# 11. Manual review template
# ============================================================
# Αυτό είναι η "τρίτη επιλογή":
# Leiden + reference evidence + DE genes → manual final decision.
#
# Εδώ θα μπορείς να συμπληρώσεις:
# - manual_decision
# - manual_final_label
# - use_for_downstream_final
#
# Δεν εφαρμόζει ακόμα labels.
# Μόνο δημιουργεί review table.

manual_template = decision_table.copy()

manual_template["manual_decision_allowed_values"] = (
    "keep_as_top_label | merge_with_same_label | subcluster_before_use | "
    "exclude_ambiguous | assign_other | custom_manual_label"
)

manual_template["manual_final_label"] = manual_template.apply(
    lambda row: row["top_reference_label"]
    if row["consensus_status"] == "high_consensus"
    else "",
    axis=1,
)

manual_template["use_for_downstream_final"] = manual_template.apply(
    lambda row: True
    if row["consensus_status"] == "high_consensus"
    and row["size_status"] == "adequate_size"
    else "",
    axis=1,
)


# ============================================================
# 12. Overview sanity table
# ============================================================

n_leiden_clusters = int(adata.obs[cluster_key].nunique())
n_reference_labels = int(adata.obs["reference_label"].nunique())

n_high = int((decision_table["consensus_status"] == "high_consensus").sum())
n_medium = int((decision_table["consensus_status"] == "medium_consensus").sum())
n_low = int((decision_table["consensus_status"] == "mixed_or_low_consensus").sum())

n_possible_split = int(
    decision_table["recommended_action"]
    .str.contains("split", case=False, na=False)
    .sum()
)

n_small = int((decision_table["size_status"] == "small_cluster_review").sum())

overview = pd.DataFrame(
    [
        {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_leiden_clusters": n_leiden_clusters,
            "n_reference_labels": n_reference_labels,
            "reference_labels_minus_leiden_clusters": n_reference_labels - n_leiden_clusters,
            "n_high_consensus_leiden_clusters": n_high,
            "n_medium_consensus_leiden_clusters": n_medium,
            "n_low_or_mixed_leiden_clusters": n_low,
            "n_possible_split_or_subcluster_clusters": n_possible_split,
            "n_small_clusters": n_small,
            "high_consensus_threshold": high_consensus_threshold,
            "medium_consensus_threshold": medium_consensus_threshold,
            "split_secondary_fraction_threshold": split_secondary_fraction_threshold,
        }
    ]
)

print("\nAnnotation sanity overview:")
print(overview.to_string(index=False))

print("\nDecision table:")
print(decision_table.to_string(index=False))

print("\nReference label spread table:")
print(reference_spread.to_string(index=False))


# ============================================================
# 13. Save tables
# ============================================================

overview_out = paths.tables / "rna_annotation_sanity_overview.csv"
decision_table_out = paths.tables / "rna_annotation_group_decision_table.csv"
reference_spread_out = paths.tables / "rna_reference_label_spread_across_leiden.csv"
top_de_out = paths.tables / "rna_annotation_top_de_genes_per_leiden.csv"
manual_template_out = paths.tables / "rna_manual_annotation_review_template.csv"

write_csv(overview, overview_out, index=False)
write_csv(decision_table, decision_table_out, index=False)
write_csv(reference_spread, reference_spread_out, index=False)
write_csv(top_de, top_de_out, index=False)
write_csv(manual_template, manual_template_out, index=False)

print("\nSaved overview to:")
print(overview_out)

print("\nSaved decision table to:")
print(decision_table_out)

print("\nSaved reference spread table to:")
print(reference_spread_out)

print("\nSaved top DE genes to:")
print(top_de_out)

print("\nSaved manual review template to:")
print(manual_template_out)


# ============================================================
# 14. Save diagnostic heatmaps
# ============================================================

figures_rna_dir = paths.figures / "rna"
figures_rna_dir.mkdir(parents=True, exist_ok=True)

heatmap_leiden_out = (
    figures_rna_dir / "pbmc_rna_sanity_leiden_reference_fraction_heatmap.png"
)

plt.figure(figsize=(14, 8))
plt.imshow(fraction_by_leiden.values, aspect="auto")
plt.xticks(
    ticks=range(fraction_by_leiden.shape[1]),
    labels=fraction_by_leiden.columns,
    rotation=90,
)
plt.yticks(
    ticks=range(fraction_by_leiden.shape[0]),
    labels=fraction_by_leiden.index,
)
plt.xlabel("Reference annotation label")
plt.ylabel("Leiden cluster")
plt.title("Reference label fractions within each Leiden cluster")
plt.colorbar(label="fraction within Leiden cluster")
plt.tight_layout()
plt.savefig(heatmap_leiden_out, dpi=300, bbox_inches="tight")
plt.close()

print("\nSaved Leiden-reference fraction heatmap to:")
print(heatmap_leiden_out)


# ============================================================
# 15. Save manifest
# ============================================================

manifest = {
    "step": "03B",
    "description": (
        "Sanity check comparing Leiden clusters, reference/model-based labels, "
        "and manual-review decision layer. Leiden clusters are treated as atomic "
        "annotation units; reference labels provide candidate biological labels."
    ),
    "inputs": {
        "reference_annotated_h5ad": str(reference_h5ad),
    },
    "outputs": {
        "strategy_counts": str(strategy_counts_out),
        "count_table": str(count_table_out),
        "fraction_by_leiden": str(fraction_by_leiden_out),
        "fraction_by_reference": str(fraction_by_reference_out),
        "overview": str(overview_out),
        "decision_table": str(decision_table_out),
        "reference_spread": str(reference_spread_out),
        "top_de": str(top_de_out),
        "manual_review_template": str(manual_template_out),
        "heatmap": str(heatmap_leiden_out),
    },
    "thresholds": {
        "high_consensus_threshold": high_consensus_threshold,
        "medium_consensus_threshold": medium_consensus_threshold,
        "split_secondary_fraction_threshold": split_secondary_fraction_threshold,
        "small_cluster_n_cells": small_cluster_n_cells,
    },
    "interpretation": {
        "default_rule": (
            "Use Leiden clusters as atomic groups. Assign a reference label only "
            "when the Leiden cluster shows sufficient consensus. Merge multiple "
            "Leiden clusters downstream if they share the same reviewed biological label."
        ),
        "do_not_do": (
            "Do not directly replace Leiden clusters with per-cell reference labels "
            "without cluster-level sanity check."
        ),
    },
}

manifest_out = (
    paths.manifests / "step03B_annotation_sanity_check_manifest.json"
)

write_json(manifest, manifest_out)

print("\nSaved manifest to:")
print(manifest_out)


# ============================================================
# 16. Τελικό summary
# ============================================================

print("\nDone.")
print("Το Step 03B annotation sanity check ολοκληρώθηκε.")
print("Άνοιξε το rna_annotation_group_decision_table.csv και το rna_manual_annotation_review_template.csv.")
print("Μετά θα φτιάξουμε Step 03C για να εφαρμόσουμε τα final reviewed groups.")