"""Βήμα 02A — Interactive RNA QC diagnostics και threshold exploration για informed preprocessing choices"""

from pathlib import Path  # ασφαλής χειρισμός paths
import sys  # για να προσθέσουμε Refactored/src στο import path
from itertools import product  # για να δοκιμάζουμε combinations από threshold lists

import pandas as pd  # για summary tables
import scanpy as sc  # για single-cell AnnData και QC metrics
import matplotlib.pyplot as plt  # για QC figures
from scipy import sparse  # για αποδοτικό χειρισμό sparse matrices
import shutil  # για backup του YAML
from datetime import datetime  # για timestamp στο backup
import yaml  # για ασφαλές γράψιμο YAML
import subprocess  # τρέχει pip από Python

subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml"])  # εγκαθιστά PyYAML

# ============================================================
# 0. Project root και imports
# ============================================================

refactor_dir = Path(__file__).resolve().parents[1]  # Refactored root folder

src_dir = refactor_dir / "src"  # shared package source folder

if str(src_dir) not in sys.path:  # αν δεν υπάρχει ήδη στο import path
    sys.path.insert(0, str(src_dir))  # επιτρέπει imports από pbmcgrn

from pbmcgrn.config import load_yaml  # YAML config loader
from pbmcgrn.io import write_csv  # reusable CSV writer
from pbmcgrn.paths import ProjectPaths  # project-aware path manager
from pbmcgrn.preprocessing.rna import calculate_rna_qc_metrics  # reusable RNA QC function


# ============================================================
# 1. Helper functions
# ============================================================

"""Βήμα 02A — Helpers για final threshold selection"""

def ask_final_single_value(prompt_text, default_value, value_type):  # ζητά χειροκίνητο single final value
    user_text = input(f"{prompt_text} final value [default: {default_value}]: ").strip()  # input από χρήστη

    if user_text == "":  # αν πατήσεις Enter
        return value_type(default_value)  # κρατά default value

    return value_type(user_text)  # επιστρέφει typed value


def choose_final_threshold_candidate(candidate_summary):  # επιλέγει ένα row από το candidate summary
    review_columns = [  # columns για καθαρή τελική επιλογή
        "candidate_id",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
        "cells_after",
        "cells_removed",
        "cell_retention_fraction",
        "genes_after",
        "genes_removed",
        "gene_retention_fraction",
    ]

    print("\n" + "=" * 80)  # separator
    print("FINAL THRESHOLD CANDIDATE SELECTION")  # section title
    print("=" * 80)  # separator

    print("\nTop candidate combinations:")  # δείχνει top combinations
    print(candidate_summary[review_columns].head(30))  # εμφανίζει top 30 rows

    print("\nΔιάλεξε ένα candidate_id από τον πίνακα.")  # οδηγία
    print("Αν θες να βάλεις τελικές τιμές χειροκίνητα, γράψε: manual")  # manual option

    selected_text = input("\nFinal candidate_id or manual [default: 1]: ").strip()  # ζητά επιλογή

    if selected_text == "":  # αν πατήσεις Enter
        selected_text = "1"  # default candidate

    if selected_text.lower() == "manual":  # manual mode
        return None  # επιστρέφει None για manual final input

    selected_id = int(selected_text)  # μετατρέπει candidate id σε int

    if selected_id not in candidate_summary["candidate_id"].values:  # έλεγχος ότι υπάρχει
        raise ValueError(f"candidate_id not found: {selected_id}")  # καθαρό error

    selected_row = candidate_summary.loc[  # βρίσκει το selected row
        candidate_summary["candidate_id"] == selected_id
    ].iloc[0]  # παίρνει τη μοναδική γραμμή

    return selected_row  # επιστρέφει το selected candidate row


"""Βήμα 02A — Helper για επιλογή ενός final candidate combination"""

def choose_final_threshold_candidate(candidate_summary):
    review_columns = [
        "candidate_id",
        "min_genes_per_cell",
        "min_cells_per_gene",
        "max_pct_mito",
        "max_total_counts",
        "cells_after",
        "cells_removed",
        "cell_retention_fraction",
        "genes_after",
        "genes_removed",
        "gene_retention_fraction",
    ]  # columns που χρειαζόμαστε για final review

    print("\n" + "=" * 80)
    print("FINAL THRESHOLD CANDIDATE SELECTION")
    print("=" * 80)

    print("\nTop candidate combinations:")
    print(candidate_summary[review_columns].head(30))  # δείχνει τα πρώτα 30 combinations

    print("\nΔιάλεξε ένα candidate_id από τον πίνακα.")
    print("Αν θες να βάλεις τελικές τιμές χειροκίνητα, γράψε: manual")

    selected_text = input("\nFinal candidate_id or manual [default: 1]: ").strip()

    if selected_text == "":
        selected_text = "1"  # default: πρώτο candidate row

    if selected_text.lower() == "manual":
        return None  # θα γίνει χειροκίνητη επιλογή μετά

    selected_id = int(selected_text)  # candidate id που διάλεξες

    if selected_id not in candidate_summary["candidate_id"].values:
        raise ValueError(f"candidate_id not found: {selected_id}")  # προστασία από λάθος id

    selected_row = (
        candidate_summary
        .loc[candidate_summary["candidate_id"] == selected_id]
        .iloc[0]
    )  # παίρνει το selected row

    return selected_row  # επιστρέφει final selected candidate


"""Βήμα 02A — Helper functions για εμφάνιση relevant statistics πριν από κάθε threshold input"""

"""Βήμα 02A — Helper για επιλογή τελικού single threshold από candidates"""

def ask_final_single_value(prompt_text, default_value, value_type):  # ζητά ένα τελικό value για YAML
    user_text = input(f"{prompt_text} final value [default: {default_value}]: ").strip()  # ζητά final value

    if user_text == "":  # αν πατήσεις Enter
        return value_type(default_value)  # κρατάει το current/default value

    return value_type(user_text)  # επιστρέφει το value στον σωστό τύπο


def show_cell_metric_statistics(adata, metric_name):  # δείχνει statistics για ένα cell-level QC metric
    stats = adata.obs[metric_name].describe()  # βασικά summary statistics

    quantiles = adata.obs[metric_name].quantile(  # useful threshold quantiles
        [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    )

    review_table = pd.DataFrame(  # compact table για decision making
        {
            "statistic": [
                "count",
                "mean",
                "std",
                "min",
                "q01",
                "q05",
                "q10",
                "q25",
                "q50",
                "q75",
                "q90",
                "q95",
                "q99",
                "max",
            ],
            "value": [
                stats["count"],
                stats["mean"],
                stats["std"],
                stats["min"],
                quantiles.loc[0.01],
                quantiles.loc[0.05],
                quantiles.loc[0.10],
                quantiles.loc[0.25],
                quantiles.loc[0.50],
                quantiles.loc[0.75],
                quantiles.loc[0.90],
                quantiles.loc[0.95],
                quantiles.loc[0.99],
                stats["max"],
            ],
        }
    )

    print("\n" + "=" * 80)  # separator
    print(f"RELEVANT STATISTICS FOR: {metric_name}")  # δείχνει ποιο metric αξιολογείς
    print("=" * 80)  # separator
    print(review_table)  # εμφανίζει τον πίνακα statistics


def show_gene_detection_statistics(adata):  # δείχνει statistics για gene-level detection
    if sparse.issparse(adata.X):  # αν το expression matrix είναι sparse
        detected_cells_per_gene = adata.X.getnnz(axis=0)  # πόσα cells εκφράζουν κάθε gene
    else:  # αν είναι dense matrix
        detected_cells_per_gene = (adata.X > 0).sum(axis=0)  # πόσα cells εκφράζουν κάθε gene

    detected_cells_per_gene = pd.Series(  # μετατροπή σε pandas Series
        detected_cells_per_gene,  # detected cell counts per gene
        index=adata.var_names,  # gene names
        name="detected_cells_per_gene",  # column name
    )

    stats = detected_cells_per_gene.describe()  # summary statistics

    quantiles = detected_cells_per_gene.quantile(  # quantiles για gene filtering
        [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
    )

    review_table = pd.DataFrame(  # compact gene-level statistics table
        {
            "statistic": [
                "count_genes",
                "mean_detected_cells",
                "std",
                "min",
                "q01",
                "q05",
                "q10",
                "q25",
                "q50",
                "q75",
                "q90",
                "q95",
                "q99",
                "max",
            ],
            "value": [
                stats["count"],
                stats["mean"],
                stats["std"],
                stats["min"],
                quantiles.loc[0.01],
                quantiles.loc[0.05],
                quantiles.loc[0.10],
                quantiles.loc[0.25],
                quantiles.loc[0.50],
                quantiles.loc[0.75],
                quantiles.loc[0.90],
                quantiles.loc[0.95],
                quantiles.loc[0.99],
                stats["max"],
            ],
        }
    )

    candidate_gene_retention = pd.DataFrame(  # δείχνει πόσα genes μένουν για common thresholds
        {
            "min_cells_per_gene": [1, 2, 3, 5, 10, 20, 50],
            "genes_retained": [
                int((detected_cells_per_gene >= threshold).sum())
                for threshold in [1, 2, 3, 5, 10, 20, 50]
            ],
        }
    )

    candidate_gene_retention["genes_removed"] = adata.n_vars - candidate_gene_retention["genes_retained"]  # removed genes
    
    candidate_gene_retention["gene_retention_fraction"] = candidate_gene_retention["genes_retained"] / adata.n_vars  # retention fraction

    print("\n" + "=" * 80)  # separator
    print("RELEVANT STATISTICS FOR: min_cells_per_gene")  # gene-level filtering section
    print("=" * 80)  # separator
    print(review_table)  # εμφανίζει statistics

    print("\nGene retention for common min_cells_per_gene thresholds:")  # candidate retention section
    print(candidate_gene_retention)  # εμφανίζει genes retained/removed


def parse_number_or_list_with_context(
    adata,
    prompt_text,
    default_value,
    value_type,
    metric_name=None,
    is_gene_level=False,
):  # εμφανίζει statistics και μετά ζητά input
    if is_gene_level:  # για min_cells_per_gene
        show_gene_detection_statistics(adata)  # δείχνει gene-level detection statistics
    elif metric_name is not None:  # για cell-level QC metrics
        show_cell_metric_statistics(adata, metric_name)  # δείχνει statistics για το συγκεκριμένο metric

    default_text = str(default_value)  # default ως text

    user_text = input(f"\n{prompt_text} [default: {default_text}]: ").strip()  # ζητά input μετά τα statistics

    if user_text == "":  # αν πατήσεις Enter
        user_text = default_text  # κρατά default

    values = [value.strip() for value in user_text.split(",")]  # επιτρέπει input τύπου 200,300,500

    parsed_values = [value_type(value) for value in values if value != ""]  # μετατρέπει values σε int/float

    return parsed_values  # επιστρέφει πάντα list
def parse_number_or_list(prompt_text, default_value, value_type):  # δέχεται έναν αριθμό ή comma-separated list
    default_text = str(default_value)  # μετατρέπει default σε string για εμφάνιση

    user_text = input(f"{prompt_text} [default: {default_text}]: ").strip()  # ζητά input από τον χρήστη

    if user_text == "":  # αν πατήσει Enter
        user_text = default_text  # χρησιμοποιεί default από YAML

    values = [value.strip() for value in user_text.split(",")]  # χωρίζει input τύπου 200,300,500

    parsed_values = [value_type(value) for value in values if value != ""]  # μετατρέπει σε int ή float

    return parsed_values  # επιστρέφει πάντα list


def count_genes_detected_in_at_least_n_cells(adata, cell_mask, min_cells_per_gene):  # μετρά genes που περνάνε min_cells
    x_subset = adata.X[cell_mask, :]  # κρατά cells που περνάνε το cell-level mask

    if sparse.issparse(x_subset):  # αν το matrix είναι sparse
        detected_cell_counts = x_subset.getnnz(axis=0)  # μετρά non-zero cells ανά gene
    else:  # αν το matrix είναι dense
        detected_cell_counts = (x_subset > 0).sum(axis=0)  # μετρά positive expression ανά gene

    n_genes_after = int((detected_cell_counts >= min_cells_per_gene).sum())  # genes που περνάνε min_cells

    return n_genes_after  # επιστρέφει αριθμό genes


def build_threshold_candidate_summary(
    adata,
    min_genes_values,
    min_cells_gene_values,
    max_pct_mito_values,
    max_total_counts_values,
):  # φτιάχνει table για όλα τα threshold combinations

    rows = []  # εδώ θα μαζέψουμε κάθε candidate result

    n_cells_before = int(adata.n_obs)  # αρχικός αριθμός cells
    n_genes_before = int(adata.n_vars)  # αρχικός αριθμός genes

    obs = adata.obs  # συντόμευση για QC metadata

    for min_genes, min_cells_gene, max_pct_mito, max_total_counts in product(  # όλα τα combinations
        min_genes_values,  # πιθανές τιμές για min genes per cell
        min_cells_gene_values,  # πιθανές τιμές για min cells per gene
        max_pct_mito_values,  # πιθανές τιμές για mito threshold
        max_total_counts_values,  # πιθανές τιμές για max counts
    ):
        min_genes_mask = obs["n_genes_by_counts"] >= min_genes  # cells που έχουν αρκετά detected genes

        mito_mask = obs["pct_counts_mt"] < max_pct_mito  # cells κάτω από mito threshold

        counts_mask = obs["total_counts"] < max_total_counts  # cells κάτω από total count threshold

        final_cell_mask = (min_genes_mask & mito_mask & counts_mask).to_numpy()  # τελικό cell-level mask

        gene_filter_cell_mask = min_genes_mask.to_numpy()  # mimic Step 02: filter genes μετά από min_genes filter

        n_cells_after = int(final_cell_mask.sum())  # τελικά cells που περνάνε thresholds

        n_genes_after = count_genes_detected_in_at_least_n_cells(  # genes που περνάνε min_cells_per_gene
            adata=adata,  # raw QC AnnData
            cell_mask=gene_filter_cell_mask,  # cells μετά το min_genes filter
            min_cells_per_gene=min_cells_gene,  # threshold για gene filtering
        )

        rows.append(  # προσθέτει candidate row
            {
                "min_genes_per_cell": int(min_genes),  # candidate min genes
                "min_cells_per_gene": int(min_cells_gene),  # candidate min cells per gene
                "max_pct_mito": float(max_pct_mito),  # candidate mito threshold
                "max_total_counts": float(max_total_counts),  # candidate max counts
                "cells_before": n_cells_before,  # αρχικά cells
                "cells_after": n_cells_after,  # cells μετά τα thresholds
                "cells_removed": n_cells_before - n_cells_after,  # removed cells
                "cell_retention_fraction": n_cells_after / n_cells_before,  # retention fraction
                "genes_before": n_genes_before,  # αρχικά genes
                "genes_after": n_genes_after,  # genes μετά το gene filter
                "genes_removed": n_genes_before - n_genes_after,  # removed genes
                "gene_retention_fraction": n_genes_after / n_genes_before,  # gene retention fraction
                "removed_by_min_genes_only": int((~min_genes_mask).sum()),  # cells που κόβονται μόνο από min genes
                "removed_by_mito_only": int((~mito_mask).sum()),  # cells που κόβονται μόνο από mito
                "removed_by_total_counts_only": int((~counts_mask).sum()),  # cells που κόβονται μόνο από total counts
            }
        )

    candidate_summary = pd.DataFrame(rows)  # μετατρέπει rows σε dataframe

    candidate_summary = candidate_summary.sort_values(  # βάζει πιο conservative/retentive combinations σε σειρά
        by=["cell_retention_fraction", "gene_retention_fraction"],  # ταξινόμηση με βάση retention
        ascending=False,  # υψηλότερο retention πρώτα
    ).reset_index(drop=True)  # καθαρό index

    return candidate_summary  # επιστρέφει candidate table

"""Βήμα 02A — QC plots με όλες τις candidate threshold γραμμές"""

def save_qc_distribution_plots(adata, outputs, threshold_candidates):  # σώζει QC plots με όλες τις τιμές που έδωσε ο χρήστης
    qc_columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]  # βασικά QC metrics

    sc.pl.violin(  # φτιάχνει violin plots για QC distributions
        adata,  # AnnData με QC metrics
        keys=qc_columns,  # metrics προς εμφάνιση
        jitter=0.4,  # δείχνει distribution των cells
        multi_panel=True,  # χωριστό panel ανά metric
        show=False,  # δεν ανοίγει interactive παράθυρο
    )

    plt.savefig(outputs["qc_violin_plot"], dpi=300, bbox_inches="tight")  # σώζει violin plot
    plt.close()  # κλείνει figure

    fig, axes = plt.subplots(1, 3, figsize=(16, 4))  # τρία histograms σε μία γραμμή

    axes[0].hist(adata.obs["n_genes_by_counts"], bins=80, color="steelblue", alpha=0.85)  # detected genes histogram
    for value in threshold_candidates["min_genes_per_cell"]:  # όλες οι candidate τιμές για min genes
        axes[0].axvline(value, linestyle="--", linewidth=1.5, label=f"{value}")  # threshold line
    axes[0].set_title("Detected genes per cell")  # τίτλος panel
    axes[0].set_xlabel("n_genes_by_counts")  # x label
    axes[0].set_ylabel("Number of cells")  # y label
    axes[0].legend(title="min genes", fontsize=7)  # legend

    axes[1].hist(adata.obs["pct_counts_mt"], bins=80, color="darkorange", alpha=0.85)  # mito histogram
    for value in threshold_candidates["max_pct_mito"]:  # όλες οι candidate τιμές για mito
        axes[1].axvline(value, linestyle="--", linewidth=1.5, label=f"{value}")  # threshold line
    axes[1].set_title("Mitochondrial percentage")  # τίτλος panel
    axes[1].set_xlabel("pct_counts_mt")  # x label
    axes[1].set_ylabel("Number of cells")  # y label
    axes[1].legend(title="max mito %", fontsize=7)  # legend

    axes[2].hist(adata.obs["total_counts"], bins=80, color="seagreen", alpha=0.85)  # total counts histogram
    for value in threshold_candidates["max_total_counts"]:  # όλες οι candidate τιμές για max counts
        axes[2].axvline(value, linestyle="--", linewidth=1.5, label=f"{value}")  # threshold line
    axes[2].set_title("Total counts per cell")  # τίτλος panel
    axes[2].set_xlabel("total_counts")  # x label
    axes[2].set_ylabel("Number of cells")  # y label
    axes[2].legend(title="max counts", fontsize=7)  # legend

    plt.tight_layout()  # καθαρό layout
    plt.savefig(outputs["qc_histogram_plot"], dpi=300, bbox_inches="tight")  # σώζει histogram plot
    plt.close()  # κλείνει figure

    sc.pl.scatter(  # scatter total counts vs mito
        adata,  # AnnData με QC metrics
        x="total_counts",  # x-axis
        y="pct_counts_mt",  # y-axis
        show=False,  # δεν εμφανίζει interactive plot
    )

    for value in threshold_candidates["max_pct_mito"]:  # όλες οι mito threshold γραμμές
        plt.axhline(value, linestyle="--", linewidth=1.2, label=f"mito {value}")  # horizontal line

    for value in threshold_candidates["max_total_counts"]:  # όλες οι count threshold γραμμές
        plt.axvline(value, linestyle=":", linewidth=1.2, label=f"counts {value}")  # vertical line

    plt.legend(fontsize=7)  # legend
    plt.savefig(outputs["qc_counts_mito_scatter"], dpi=300, bbox_inches="tight")  # σώζει scatter
    plt.close()  # κλείνει figure

    sc.pl.scatter(  # scatter total counts vs detected genes
        adata,  # AnnData με QC metrics
        x="total_counts",  # x-axis
        y="n_genes_by_counts",  # y-axis
        show=False,  # δεν εμφανίζει interactive plot
    )

    for value in threshold_candidates["min_genes_per_cell"]:  # όλες οι min genes γραμμές
        plt.axhline(value, linestyle="--", linewidth=1.2, label=f"genes {value}")  # horizontal line

    for value in threshold_candidates["max_total_counts"]:  # όλες οι count threshold γραμμές
        plt.axvline(value, linestyle=":", linewidth=1.2, label=f"counts {value}")  # vertical line

    plt.legend(fontsize=7)  # legend
    plt.savefig(outputs["qc_counts_genes_scatter"], dpi=300, bbox_inches="tight")  # σώζει scatter
    plt.close()  # κλείνει figure

# ============================================================
# 2. Initialize project-aware paths
# ============================================================

paths = ProjectPaths(refactor_root=refactor_dir)  # διαβάζει active project από registry
paths.ensure_output_dirs()  # δημιουργεί output dirs αν λείπουν

print("\nRefactored project root:")  # ενημερωτικό print
print(paths.refactor_root)  # δείχνει Refactored root

print("\nActive project id:")  # ενημερωτικό print
print(paths.project_id)  # δείχνει active project id

print("\nActive project root:")  # ενημερωτικό print
print(paths.project_root)  # δείχνει active project root


# ============================================================
# 3. Load config
# ============================================================

config_path = paths.configs / "preprocessing_rna.yaml"  # ίδιο YAML με Step 02

print("\nLoading RNA preprocessing config:")  # ενημερωτικό print
print(config_path)  # δείχνει YAML path

config = load_yaml(config_path)  # φορτώνει YAML

rna_cfg = config["rna_preprocessing"]  # QC/filtering config
qc_diag_cfg = config["qc_diagnostics"]  # QC diagnostics config
qc_diag_outputs_cfg = qc_diag_cfg["outputs"]  # diagnostics output paths

outputs = {  # project-aware resolved output paths
    key: paths.prepare_output_path(value, label=key)  # μετατρέπει relative project path σε absolute path
    for key, value in qc_diag_outputs_cfg.items()  # διαβάζει όλα τα diagnostics outputs από YAML
}

print("\nResolved QC diagnostics outputs:")  # ενημερωτικό print
for key, value in outputs.items():  # loop στα output paths
    print(f"{key}: {value}")  # δείχνει κάθε output path


# ============================================================
# 4. Load raw RNA and calculate QC metrics
# ============================================================

rna_raw_path = paths.processed_rna / "pbmc_rna_raw.h5ad"  # raw RNA object από Step 01

print("\nLoading raw RNA object:")  # ενημερωτικό print
print(rna_raw_path)  # δείχνει input path

if not rna_raw_path.exists():  # αν λείπει το raw RNA
    raise FileNotFoundError(f"Raw RNA h5ad not found:\n{rna_raw_path}")  # σταματά με καθαρό error

rna_raw = sc.read_h5ad(rna_raw_path)  # φορτώνει raw RNA AnnData

print("\nRaw RNA object:")  # ενημερωτικό print
print(rna_raw)  # δείχνει n_cells × n_genes

default_mito_prefix = rna_cfg["mito_prefix"]  # default mito prefix από YAML

mito_prefix_input = input(f"\nmito_prefix [default: {default_mito_prefix}]: ").strip()  # ζητά prefix

mito_prefix = mito_prefix_input if mito_prefix_input else default_mito_prefix  # κρατά default αν πατήσεις Enter

rna_qc = calculate_rna_qc_metrics(  # υπολογίζει QC metrics
    adata=rna_raw,  # raw RNA AnnData
    mito_prefix=mito_prefix,  # prefix που επέλεξες
)

n_mito_genes = int(rna_qc.var["mt"].sum())  # μετρά mitochondrial genes που βρέθηκαν

print("\nMitochondrial genes detected:")  # ενημερωτικό print
print(f"{n_mito_genes} genes with prefix '{mito_prefix}'")  # δείχνει αν το prefix είναι σωστό


# ============================================================
# 5. Save raw QC statistics
# ============================================================

qc_columns = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]  # QC metrics για statistics

qc_summary = rna_qc.obs[qc_columns].describe().T.reset_index(names="metric")  # summary statistics

qc_quantiles = rna_qc.obs[qc_columns].quantile(qc_diag_cfg["quantiles"]).T  # quantile table

qc_quantiles = qc_quantiles.reset_index(names="metric")  # κάνει το metric κανονική στήλη

write_csv(qc_summary, outputs["raw_qc_summary"], index=False)  # σώζει raw QC summary

write_csv(qc_quantiles, outputs["raw_qc_quantiles"], index=False)  # σώζει raw QC quantiles

print("\nRaw QC summary:")  # ενημερωτικό print
print(qc_summary)  # δείχνει summary statistics

print("\nRaw QC quantiles:")  # ενημερωτικό print
print(qc_quantiles)  # δείχνει quantile table





# ============================================================
# 7. Ask user for threshold candidates
# ============================================================

print("\nNow enter threshold candidates.")  # ενημερωτικό print
print("Μπορείς να δώσεις έναν αριθμό ή λίστα με κόμμα, π.χ. 200 ή 200,300,500.")  # οδηγία input

"""Βήμα 02A — Threshold inputs με relevant statistics ακριβώς πριν από κάθε απόφαση"""

min_genes_values = parse_number_or_list_with_context(
    adata=rna_qc,  # raw AnnData με QC metrics
    prompt_text="min_genes_per_cell",  # threshold για χαμηλά detected genes
    default_value=rna_cfg["min_genes_per_cell"],  # default από YAML
    value_type=int,  # integer input
    metric_name="n_genes_by_counts",  # relevant statistic για αυτό το threshold
)  # δείχνει n_genes_by_counts statistics και μετά ζητά input

min_cells_gene_values = parse_number_or_list_with_context(
    adata=rna_qc,  # raw AnnData με QC metrics
    prompt_text="min_cells_per_gene",  # threshold για gene-level detection
    default_value=rna_cfg["min_cells_per_gene"],  # default από YAML
    value_type=int,  # integer input
    is_gene_level=True,  # ειδικό gene-level statistics table
)  # δείχνει detected-cells-per-gene statistics και μετά ζητά input

max_pct_mito_values = parse_number_or_list_with_context(
    adata=rna_qc,  # raw AnnData με QC metrics
    prompt_text="max_pct_mito",  # threshold για high mitochondrial percentage
    default_value=rna_cfg["max_pct_mito"],  # default από YAML
    value_type=float,  # float input
    metric_name="pct_counts_mt",  # relevant statistic για mitochondrial filtering
)  # δείχνει pct_counts_mt statistics και μετά ζητά input

max_total_counts_values = parse_number_or_list_with_context(
    adata=rna_qc,  # raw AnnData με QC metrics
    prompt_text="max_total_counts",  # threshold για extreme library size / possible doublets
    default_value=rna_cfg["max_total_counts"],  # default από YAML
    value_type=float,  # float input
    metric_name="total_counts",  # relevant statistic για total-count filtering
)  # δείχνει total_counts statistics και μετά ζητά input

threshold_candidates = {  # κρατά όλες τις candidate τιμές που έδωσε ο χρήστης
    "min_genes_per_cell": min_genes_values,  # λίστα candidate min genes
    "min_cells_per_gene": min_cells_gene_values,  # λίστα candidate min cells per gene
    "max_pct_mito": max_pct_mito_values,  # λίστα candidate mito thresholds
    "max_total_counts": max_total_counts_values,  # λίστα candidate max counts
}

save_qc_distribution_plots(  # σώζει plots με όλες τις candidate γραμμές
    adata=rna_qc,  # raw AnnData με QC metrics
    outputs=outputs,  # resolved output paths
    threshold_candidates=threshold_candidates,  # candidate thresholds από user input
)

print("\nSaved QC distribution plots with user-provided threshold candidates.")  # επιβεβαίωση
# ============================================================
# 8. Evaluate threshold candidates
# ============================================================

"""Βήμα 02A — Build threshold candidate summary και πρόσθεσε candidate_id"""

candidate_summary = build_threshold_candidate_summary(
    adata=rna_qc,  # raw AnnData με QC metrics
    min_genes_values=min_genes_values,  # min genes list
    min_cells_gene_values=min_cells_gene_values,  # min cells per gene list
    max_pct_mito_values=max_pct_mito_values,  # mito list
    max_total_counts_values=max_total_counts_values,  # max counts list
)  # φτιάχνει candidate table για όλα τα threshold combinations

if "candidate_id" not in candidate_summary.columns:  # αποφεύγει duplicate candidate_id
    candidate_summary.insert(
        0,  # πρώτη στήλη
        "candidate_id",  # όνομα στήλης
        range(1, candidate_summary.shape[0] + 1),  # ids 1, 2, 3...
    )  # δίνει μοναδικό id σε κάθε threshold combination

write_csv(
    candidate_summary,  # candidate table με candidate_id
    outputs["threshold_candidate_summary"],  # output path από YAML
    index=False,  # δεν σώζει pandas index
)  # σώζει candidate summary

print("\nThreshold candidate summary saved to:")  # ενημερωτικό print
print(outputs["threshold_candidate_summary"])  # δείχνει output path

print("\nTop candidate combinations by cell/gene retention:")  # ενημερωτικό print
print(candidate_summary.head(20))  # δείχνει τα πρώτα 20 combinations

print("\nTop candidate combinations by cell/gene retention:")  # ενημερωτικό print
print(candidate_summary.head(20))  # δείχνει τα πρώτα 20 combinations
"""Βήμα 02A — Εμφάνιση QC summary statistics, candidate thresholds και saved plots"""

pd.set_option("display.max_columns", None)  # δείχνει όλες τις στήλες στα printed pandas tables
pd.set_option("display.width", 200)  # κάνει πιο πλατιά την εκτύπωση στο console
pd.set_option("display.max_rows", 30)  # δείχνει αρκετές γραμμές χωρίς να γεμίζει υπερβολικά το console


print("\n" + "=" * 80)  # separator για πιο καθαρό output
print("RAW RNA QC SUMMARY STATISTICS")  # τίτλος section
print("=" * 80)  # separator

print(qc_summary)  # δείχνει count, mean, std, min, quartiles, max για κάθε QC metric


print("\n" + "=" * 80)  # separator
print("RAW RNA QC QUANTILES")  # τίτλος section
print("=" * 80)  # separator

print(qc_quantiles)  # δείχνει τα quantiles που βοηθούν στην επιλογή thresholds


print("\n" + "=" * 80)  # separator
print("TOP THRESHOLD CANDIDATE COMBINATIONS")  # τίτλος section
print("=" * 80)  # separator

columns_to_review = [  # κρατάμε τις πιο χρήσιμες στήλες για γρήγορη αξιολόγηση
    "min_genes_per_cell",  # candidate threshold για detected genes per cell
    "min_cells_per_gene",  # candidate threshold για detected cells per gene
    "max_pct_mito",  # candidate mitochondrial threshold
    "max_total_counts",  # candidate total count threshold
    "cells_after",  # πόσα cells μένουν
    "cells_removed",  # πόσα cells αφαιρούνται
    "cell_retention_fraction",  # ποσοστό cells που μένει
    "genes_after",  # πόσα genes μένουν
    "genes_removed",  # πόσα genes αφαιρούνται
    "gene_retention_fraction",  # ποσοστό genes που μένει
]

print(candidate_summary[columns_to_review].head(25))  # δείχνει τα top 25 combinations


print("\n" + "=" * 80)  # separator
print("CURRENT YAML THRESHOLDS")  # τίτλος section
print("=" * 80)  # separator

print(f"mito_prefix: {mito_prefix}")  # mitochondrial prefix που χρησιμοποιήθηκε
print(f"min_genes_per_cell: {rna_cfg['min_genes_per_cell']}")  # current min genes από YAML
print(f"min_cells_per_gene: {rna_cfg['min_cells_per_gene']}")  # current min cells per gene από YAML
print(f"max_pct_mito: {rna_cfg['max_pct_mito']}")  # current mito threshold από YAML
print(f"max_total_counts: {rna_cfg['max_total_counts']}")  # current max counts από YAML


print("\n" + "=" * 80)  # separator
print("OPENING SAVED QC PLOTS")  # τίτλος section
print("=" * 80)  # separator

plot_paths = {  # plots που αποθηκεύτηκαν στο Step 02A
    "QC violin plot": outputs["qc_violin_plot"],  # violin distributions
    "QC histogram plot": outputs["qc_histogram_plot"],  # histograms με threshold lines
    "Total counts vs pct mito": outputs["qc_counts_mito_scatter"],  # scatter counts/mito
    "Total counts vs detected genes": outputs["qc_counts_genes_scatter"],  # scatter counts/genes
}

for plot_title, plot_path in plot_paths.items():  # loop σε κάθε saved plot
    image = plt.imread(plot_path)  # διαβάζει το saved PNG

    plt.figure(figsize=(12, 7))  # μέγεθος figure για εμφάνιση
    plt.imshow(image)  # εμφανίζει την εικόνα
    plt.axis("off")  # κρύβει axes γιατί είναι ήδη plot image
    plt.title(plot_title)  # βάζει τίτλο στο displayed figure
    plt.show()  # εμφανίζει το plot στο Spyder/Jupyter

    print(f"{plot_title}: {plot_path}")  # τυπώνει και το path για να το βρεις στον φάκελο reports


print("\n" + "=" * 80)  # separator
print("HOW TO USE THESE RESULTS")  # μικρή οδηγία
print("=" * 80)  # separator

print("1. Κοίτα πρώτα τα raw QC histograms και quantiles.")  # οδηγία για informed choice
print("2. Δες πόσα cells/genes κρατά κάθε candidate threshold combination.")  # οδηγία για candidate table
print("3. Διάλεξε ένα single final value για κάθε threshold.")  # τελικό YAML πρέπει να έχει μία τιμή ανά parameter
print("4. Βάλε τα selected values στο preprocessing_rna.yaml και ξανατρέξε Step 02.")  # reproducible final run

"""Βήμα 02A — Optional update του preprocessing_rna.yaml με τα τελικά selected RNA QC thresholds"""

print("\n" + "=" * 80)  # separator
print("OPTIONAL: UPDATE preprocessing_rna.yaml WITH FINAL SELECTED QC THRESHOLDS")  # section title
print("=" * 80)  # separator

print("Current YAML values:")  # δείχνει τις τωρινές τιμές
print(f"mito_prefix: {rna_cfg['mito_prefix']}")  # current mito prefix
print(f"min_genes_per_cell: {rna_cfg['min_genes_per_cell']}")  # current min genes
print(f"min_cells_per_gene: {rna_cfg['min_cells_per_gene']}")  # current min cells per gene
print(f"max_pct_mito: {rna_cfg['max_pct_mito']}")  # current max mito
print(f"max_total_counts: {rna_cfg['max_total_counts']}")  # current max total counts

apply_update = input("\nUpdate preprocessing_rna.yaml with final selected values? [y/N]: ").strip().lower()  # safety prompt

"""Βήμα 02A — Update preprocessing_rna.yaml με selected candidate ή manual final values"""

if apply_update == "y":  # ενημερώνει YAML μόνο αν το ζητήσεις
    selected_candidate = choose_final_threshold_candidate(candidate_summary)  # επιλέγεις ένα candidate row ή manual

    final_mito_prefix = input(f"mito_prefix final value [default: {mito_prefix}]: ").strip()  # final mito prefix

    if final_mito_prefix == "":
        final_mito_prefix = mito_prefix  # κρατά το prefix που χρησιμοποιήθηκε στο 02A

    if selected_candidate is not None:
        final_min_genes_per_cell = int(selected_candidate["min_genes_per_cell"])  # από selected row
        final_min_cells_per_gene = int(selected_candidate["min_cells_per_gene"])  # από selected row
        final_max_pct_mito = float(selected_candidate["max_pct_mito"])  # από selected row
        final_max_total_counts = float(selected_candidate["max_total_counts"])  # από selected row

        print("\nSelected final threshold combination:")
        print(selected_candidate[
            [
                "candidate_id",
                "min_genes_per_cell",
                "min_cells_per_gene",
                "max_pct_mito",
                "max_total_counts",
                "cells_after",
                "cells_removed",
                "cell_retention_fraction",
                "genes_after",
                "genes_removed",
                "gene_retention_fraction",
            ]
        ])

    else:
        final_min_genes_per_cell = ask_final_single_value(
            prompt_text="min_genes_per_cell",
            default_value=rna_cfg["min_genes_per_cell"],
            value_type=int,
        )  # manual final min genes

        final_min_cells_per_gene = ask_final_single_value(
            prompt_text="min_cells_per_gene",
            default_value=rna_cfg["min_cells_per_gene"],
            value_type=int,
        )  # manual final min cells per gene

        final_max_pct_mito = ask_final_single_value(
            prompt_text="max_pct_mito",
            default_value=rna_cfg["max_pct_mito"],
            value_type=float,
        )  # manual final mito threshold

        final_max_total_counts = ask_final_single_value(
            prompt_text="max_total_counts",
            default_value=rna_cfg["max_total_counts"],
            value_type=float,
        )  # manual final max counts

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # timestamp για backup
    backup_config_path = config_path.with_suffix(f".backup_{timestamp}.yaml")  # backup YAML path

    shutil.copy2(config_path, backup_config_path)  # backup πριν το overwrite

    config["rna_preprocessing"]["mito_prefix"] = final_mito_prefix  # final mito prefix
    config["rna_preprocessing"]["min_genes_per_cell"] = final_min_genes_per_cell  # final min genes
    config["rna_preprocessing"]["min_cells_per_gene"] = final_min_cells_per_gene  # final min cells per gene
    config["rna_preprocessing"]["max_pct_mito"] = final_max_pct_mito  # final max mito
    config["rna_preprocessing"]["max_total_counts"] = final_max_total_counts  # final max counts

    with open(config_path, "w", encoding="utf-8") as yaml_file:
        yaml.safe_dump(
            config,
            yaml_file,
            sort_keys=False,
            allow_unicode=True,
        )  # γράφει updated YAML

    selected_thresholds = pd.DataFrame(
        [
            {
                "mito_prefix": final_mito_prefix,
                "min_genes_per_cell": final_min_genes_per_cell,
                "min_cells_per_gene": final_min_cells_per_gene,
                "max_pct_mito": final_max_pct_mito,
                "max_total_counts": final_max_total_counts,
                "source_script": "run_02A_rna_qc_diagnostics.py",
                "yaml_updated": str(config_path),
                "yaml_backup": str(backup_config_path),
            }
        ]
    )  # record των τελικών selected thresholds

    selected_thresholds_out = paths.prepare_output_path(
        "reports/tables/rna/rna_selected_qc_thresholds.csv",
        label="rna_selected_qc_thresholds",
    )  # output table

    write_csv(selected_thresholds, selected_thresholds_out, index=False)  # σώζει final thresholds

    print("\nUpdated preprocessing_rna.yaml:")
    print(config_path)

    print("\nBackup saved to:")
    print(backup_config_path)

    print("\nSelected thresholds saved to:")
    print(selected_thresholds_out)

    print("\nFinal selected values written to YAML:")
    print(selected_thresholds)

else:
    print("\nSkipped YAML update.")
    print("Μπορείς να ξανατρέξεις το 02A ή να κάνεις manual edit στο preprocessing_rna.yaml.")  # guidance
# ============================================================
# 9. Final guidance
# ============================================================

print("\nDone.")  # τελικό print
print("Το Step 02A RNA QC diagnostics ολοκληρώθηκε.")  # επιβεβαίωση ολοκλήρωσης
print("Άνοιξε τα QC plots και το threshold_candidate_summary.csv.")  # οδηγία επόμενου βήματος
print("Μετά διάλεξε ένα τελικό single value για κάθε threshold και βάλε το στο preprocessing_rna.yaml.")  # reproducibility guidance

"""Βήμα 02A — Optional YAML update με ένα selected final threshold combination"""

print("\n" + "=" * 80)  # separator
print("OPTIONAL: UPDATE preprocessing_rna.yaml WITH FINAL SELECTED QC THRESHOLDS")  # section title
print("=" * 80)  # separator

print("Current YAML values:")  # εμφανίζει τωρινές YAML τιμές
print(f"mito_prefix: {rna_cfg['mito_prefix']}")  # current mito prefix
print(f"min_genes_per_cell: {rna_cfg['min_genes_per_cell']}")  # current min genes
print(f"min_cells_per_gene: {rna_cfg['min_cells_per_gene']}")  # current min cells per gene
print(f"max_pct_mito: {rna_cfg['max_pct_mito']}")  # current mito threshold
print(f"max_total_counts: {rna_cfg['max_total_counts']}")  # current max counts

apply_update = input("\nUpdate preprocessing_rna.yaml with final selected values? [y/N]: ").strip().lower()  # confirmation

if apply_update == "y":  # συνεχίζει μόνο αν γράψεις y
    selected_candidate = choose_final_threshold_candidate(candidate_summary)  # επιλέγει candidate row ή manual mode

    final_mito_prefix = input(f"mito_prefix final value [default: {mito_prefix}]: ").strip()  # final mito prefix

    if final_mito_prefix == "":  # αν πατήσεις Enter
        final_mito_prefix = mito_prefix  # κρατά το prefix που χρησιμοποιήθηκε στο diagnostics

    if selected_candidate is not None:  # αν διάλεξες candidate_id
        final_min_genes_per_cell = int(selected_candidate["min_genes_per_cell"])  # final min genes από row
        final_min_cells_per_gene = int(selected_candidate["min_cells_per_gene"])  # final min cells per gene από row
        final_max_pct_mito = float(selected_candidate["max_pct_mito"])  # final mito από row
        final_max_total_counts = float(selected_candidate["max_total_counts"])  # final max counts από row

        print("\nSelected final threshold combination:")  # confirmation
        print(
            selected_candidate[
                [
                    "candidate_id",
                    "min_genes_per_cell",
                    "min_cells_per_gene",
                    "max_pct_mito",
                    "max_total_counts",
                    "cells_after",
                    "cells_removed",
                    "cell_retention_fraction",
                    "genes_after",
                    "genes_removed",
                    "gene_retention_fraction",
                ]
            ]
        )  # δείχνει το selected row

    else:  # manual mode
        final_min_genes_per_cell = ask_final_single_value(
            prompt_text="min_genes_per_cell",  # parameter name
            default_value=rna_cfg["min_genes_per_cell"],  # current YAML default
            value_type=int,  # integer
        )  # manual final min genes

        final_min_cells_per_gene = ask_final_single_value(
            prompt_text="min_cells_per_gene",  # parameter name
            default_value=rna_cfg["min_cells_per_gene"],  # current YAML default
            value_type=int,  # integer
        )  # manual final min cells per gene

        final_max_pct_mito = ask_final_single_value(
            prompt_text="max_pct_mito",  # parameter name
            default_value=rna_cfg["max_pct_mito"],  # current YAML default
            value_type=float,  # float
        )  # manual final mito threshold

        final_max_total_counts = ask_final_single_value(
            prompt_text="max_total_counts",  # parameter name
            default_value=rna_cfg["max_total_counts"],  # current YAML default
            value_type=float,  # float
        )  # manual final max counts

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # timestamp για backup
    backup_config_path = config_path.with_suffix(f".backup_{timestamp}.yaml")  # backup filename

    shutil.copy2(config_path, backup_config_path)  # κάνει backup πριν το overwrite

    config["rna_preprocessing"]["mito_prefix"] = final_mito_prefix  # γράφει final mito prefix
    config["rna_preprocessing"]["min_genes_per_cell"] = final_min_genes_per_cell  # γράφει final min genes
    config["rna_preprocessing"]["min_cells_per_gene"] = final_min_cells_per_gene  # γράφει final min cells/gene
    config["rna_preprocessing"]["max_pct_mito"] = final_max_pct_mito  # γράφει final max mito
    config["rna_preprocessing"]["max_total_counts"] = final_max_total_counts  # γράφει final max counts

    with open(config_path, "w", encoding="utf-8") as yaml_file:  # ανοίγει YAML για overwrite
        yaml.safe_dump(
            config,  # updated config
            yaml_file,  # file handle
            sort_keys=False,  # κρατά τη σειρά όσο γίνεται
            allow_unicode=True,  # επιτρέπει ελληνικά
        )  # γράφει updated YAML

    selected_thresholds = pd.DataFrame(
        [
            {
                "mito_prefix": final_mito_prefix,  # final mito prefix
                "min_genes_per_cell": final_min_genes_per_cell,  # final min genes
                "min_cells_per_gene": final_min_cells_per_gene,  # final min cells per gene
                "max_pct_mito": final_max_pct_mito,  # final max mito
                "max_total_counts": final_max_total_counts,  # final max counts
                "source_script": "run_02A_rna_qc_diagnostics.py",  # source script
                "yaml_updated": str(config_path),  # updated YAML path
                "yaml_backup": str(backup_config_path),  # backup YAML path
            }
        ]
    )  # table με final selected thresholds

    selected_thresholds_out = paths.prepare_output_path(
        "reports/tables/rna/rna_selected_qc_thresholds.csv",  # project-local output
        label="rna_selected_qc_thresholds",  # path label
    )  # resolved output path

    write_csv(selected_thresholds, selected_thresholds_out, index=False)  # σώζει selected thresholds table

    print("\nUpdated preprocessing_rna.yaml:")  # confirmation
    print(config_path)  # updated YAML path

    print("\nBackup saved to:")  # confirmation
    print(backup_config_path)  # backup path

    print("\nSelected thresholds saved to:")  # confirmation
    print(selected_thresholds_out)  # selected thresholds csv path

    print("\nFinal selected values written to YAML:")  # confirmation
    print(selected_thresholds)  # δείχνει final selected values

else:  # αν δεν γράψεις y
    print("\nSkipped YAML update.")  # confirmation
    print("Μπορείς να ξανατρέξεις το 02A ή να κάνεις manual edit στο preprocessing_rna.yaml.")  # guidance