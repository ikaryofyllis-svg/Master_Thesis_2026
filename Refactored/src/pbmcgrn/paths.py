"""Βήμα 00Q — Project-aware path management για καθαρά outputs ανά project"""

from dataclasses import dataclass  # compact class για centralized path management
from pathlib import Path  # ασφαλής χειρισμός αρχείων και φακέλων
from typing import Any, Dict, Optional  # type hints για καθαρότερο code
import yaml  # YAML loading για project_registry.yaml


# ============================================================
# 1. Helper function για loading YAML
# ============================================================

def _load_yaml(path):
    """
    Small internal YAML loader.

    Avoids circular imports from pbmcgrn.config.
    """

    path = Path(path)  # μετατροπή input σε Path object

    if not path.exists():
        raise FileNotFoundError(f"YAML file not found: {path}")  # καθαρό error αν λείπει το YAML

    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)  # διαβάζει YAML σε Python dict

    return data or {}  # αν το YAML είναι άδειο, επιστρέφει άδειο dict


# ============================================================
# 2. Project-aware path class
# ============================================================

@dataclass
class ProjectPaths:
    """
    Centralized path management for the refactored GRN project.

    Core idea:
    - refactor_root contains shared code and scripts.
    - each dataset/project has its own workspace under:
        Refactored/projects/<project_id>/

    This prevents mixing tables, figures, h5ad files, models and manifests
    from different datasets.
    """

    refactor_root: Path  # Refactored project root
    project_id: Optional[str] = None  # optional override, αλλιώς διαβάζεται από project_registry.yaml

    def __post_init__(self):
        """
        Initialize all project-aware paths.
        """

        self.refactor_root = Path(self.refactor_root).resolve()  # absolute path για Refactored root

        self.legacy_root = self.refactor_root.parent  # αρχικό MSc thesis root, read-only reference
        self.project_registry = self.refactor_root / "project_registry.yaml"  # global project registry

        registry = self._load_project_registry()  # φορτώνει registry ή fallback

        if self.project_id is None:
            self.project_id = self._get_active_project_id(registry)  # παίρνει active project id από YAML

        self.project_metadata = self._get_project_metadata(
            registry=registry,
            project_id=self.project_id,
        )  # metadata για το active project

        project_root_string = self.project_metadata.get(
            "project_root",
            f"projects/{self.project_id}",
        )  # fallback αν λείπει project_root από registry

        self.project_root = self.resolve_refactor_path(project_root_string)  # λύνει project_root από Refactored

        self.projects_root = self.refactor_root / "projects"  # default parent folder για projects

        self.src = self.refactor_root / "src"  # shared Python package code
        self.scripts = self.refactor_root / "scripts"  # shared scripts

        self.configs = self.project_root / "configs"  # project-local config files

        self.data = self.project_root / "data"  # project-local data root
        self.raw = self.data / "raw"  # raw files copied/linked for this project
        self.external = self.data / "external"  # external resources
        self.interim = self.data / "interim"  # intermediate outputs
        self.processed = self.data / "processed"  # processed project outputs

        self.processed_rna = self.processed / "rna"  # processed RNA h5ad files
        self.processed_atac = self.processed / "atac"  # processed ATAC files
        self.processed_multiome = self.processed / "multiome"
        self.processed_integration = self.processed / "integration"  # pseudo-bulk integration outputs
        self.processed_grn = self.processed / "grn"  # candidate GRN outputs
        self.processed_gnn = self.processed / "gnn"  # GNN-ready outputs

        self.reports = self.project_root / "reports"  # project-local reports root

        self.tables = self.reports / "tables"  # all report tables
        self.tables_rna = self.tables / "rna"  # RNA-specific tables
        self.tables_atac = self.tables / "atac"  # ATAC-specific tables
        self.tables_annotation = self.tables / "annotation"  # annotation sanity/review tables
        self.tables_integration = self.tables / "integration"  # integration tables
        self.tables_grn = self.tables / "grn"  # GRN tables
        self.tables_gnn = self.tables / "gnn"  # GNN tables

        self.figures = self.reports / "figures"  # all figures
        self.figures_rna = self.figures / "rna"  # RNA figures
        self.figures_atac = self.figures / "atac"  # ATAC figures
        self.figures_annotation = self.figures / "annotation"  # annotation figures
        self.figures_integration = self.figures / "integration"  # integration figures
        self.figures_grn = self.figures / "grn"  # GRN figures
        self.figures_gnn = self.figures / "gnn"  # GNN figures

        self.manifests = self.reports / "manifests"  # reproducibility manifests
        self.logs = self.reports / "logs"  # logs

        self.models = self.project_root / "models"  # project model artifacts
        self.models_annotation = self.models / "annotation"  # annotation model artifacts
        self.models_gnn = self.models / "gnn"  # GNN model artifacts

        self.notebooks = self.project_root / "notebooks"  # project notebooks

        self.annotation_tables = self.tables_annotation  # alias για readable Step 03B code
        self.annotation_figures = self.figures_annotation  # alias για readable Step 03B code

        self.legacy_data = self.legacy_root / "data"  # old data folder, read-only by default
        self.legacy_raw = self.legacy_data / "raw"  # old raw folder
        self.legacy_processed = self.legacy_data / "processed"  # old processed folder
        self.legacy_results = self.legacy_root / "results"  # old results folder

    # ========================================================
    # Registry helpers
    # ========================================================

    def _load_project_registry(self) -> Dict[str, Any]:
        """
        Load Refactored/project_registry.yaml.

        If missing, return a minimal fallback registry for pbmc_10k_v1.
        """

        if not self.project_registry.exists():
            return {
                "active_project_id": "pbmc_10k_v1",
                "projects": {
                    "pbmc_10k_v1": {
                        "project_root": "projects/pbmc_10k_v1",
                    }
                },
            }  # fallback για backward compatibility

        return _load_yaml(self.project_registry)  # κανονική φόρτωση project registry

    def _get_active_project_id(self, registry: Optional[Dict[str, Any]] = None) -> str:
        """
        Read active project id from loaded registry.
        """

        registry = registry or self._load_project_registry()  # reuse loaded registry αν υπάρχει

        if "active_project_id" not in registry:
            raise KeyError("project_registry.yaml exists but has no active_project_id")  # καθαρό error

        return registry["active_project_id"]  # π.χ. pbmc_10k_v1

    def _get_project_metadata(self, registry: Dict[str, Any], project_id: str) -> Dict[str, Any]:
        """
        Return metadata for the requested project id.
        """

        projects = registry.get("projects", {})  # παίρνει όλα τα project records

        if project_id not in projects:
            raise KeyError(
                f"Project id '{project_id}' was not found in project_registry.yaml"
            )  # σταματάει αν το project δεν υπάρχει στο registry

        return projects[project_id] or {}  # επιστρέφει metadata dict

    # ========================================================
    # Directory creation
    # ========================================================

    def ensure_output_dirs(self):
        """
        Create all standard project output folders.
        """

        folders = [
            self.project_root,

            self.configs,

            self.raw,
            self.external,
            self.interim,
            self.processed,
            self.processed_rna,
            self.processed_atac,
            self.processed_multiome,
            self.processed_integration,
            self.processed_grn,
            self.processed_gnn,

            self.reports,
            self.tables,
            self.tables_rna,
            self.tables_atac,
            self.tables_annotation,
            self.tables_integration,
            self.tables_grn,
            self.tables_gnn,

            self.figures,
            self.figures_rna,
            self.figures_atac,
            self.figures_annotation,
            self.figures_integration,
            self.figures_grn,
            self.figures_gnn,

            self.manifests,
            self.logs,

            self.models,
            self.models_annotation,
            self.models_gnn,

            self.notebooks,
        ]  # όλοι οι standard project folders

        for folder in folders:
            folder.mkdir(parents=True, exist_ok=True)  # δημιουργεί folder αν λείπει

    # ========================================================
    # Resolve paths from config
    # ========================================================

    def resolve_project_path(self, path_string):
        """
        Resolve a path written in a project config.

        Rules:
        - Absolute path: return as-is.
        - Relative path: resolve from active project root.

        Example:
        data/processed/rna/file.h5ad
        -> Refactored/projects/<project_id>/data/processed/rna/file.h5ad
        """

        if path_string is None:
            raise ValueError("Cannot resolve None as a project path")  # προστασία από άδειο config value

        path = Path(path_string)  # μετατροπή config value σε Path

        if path.is_absolute():
            return path.resolve()  # absolute paths μένουν absolute

        return (self.project_root / path).resolve()  # relative paths λύνονται από active project root


    def resolve_refactor_path(self, path_string):
        """
        Resolve path from Refactored root.

        Use this only for shared code/global files.
        Project outputs should use resolve_project_path.
        """

        if path_string is None:
            raise ValueError("Cannot resolve None as a refactor path")  # προστασία από άδειο config value

        path = Path(path_string)  # μετατροπή config value σε Path

        if path.is_absolute():
            return path.resolve()  # absolute paths μένουν absolute

        return (self.refactor_root / path).resolve()  # relative paths λύνονται από Refactored root

    # ========================================================
    # Output safety guards
    # ========================================================

    def assert_inside_project(self, path, label="path"):
        """
        Ensure an output path stays inside the active project workspace.

        Use this for outputs, not necessarily for raw external inputs.
        """

        path = Path(path).resolve()  # absolute candidate path
        project_root = Path(self.project_root).resolve()  # absolute active project root

        try:
            path.relative_to(project_root)  # πετυχαίνει μόνο αν path είναι μέσα στο project
        except ValueError:
            raise ValueError(
                f"{label} is outside the active project root.\n"
                f"Path: {path}\n"
                f"Project root: {project_root}\n"
                "Use paths.resolve_project_path(...) for project outputs."
            )  # σταματάει για να μη γράψεις σε global Refactored/data ή Refactored/reports

        return path  # επιστρέφει validated path

    def prepare_output_path(self, path_string, label="output"):
        """
        Resolve, validate and create parent directory for an output path.

        This is the safest helper for project outputs from YAML configs.
        """

        path = self.resolve_project_path(path_string)  # λύνει relative output από active project root
        path = self.assert_inside_project(path, label=label)  # ελέγχει ότι μένει μέσα στο project
        path.parent.mkdir(parents=True, exist_ok=True)  # δημιουργεί parent folder
        return path  # επιστρέφει έτοιμο output path

    # ========================================================
    # Summary
    # ========================================================

    def as_dict(self):
        """
        Export key paths as dictionary for manifests/debugging.
        """

        return {
            "refactor_root": str(self.refactor_root),
            "project_id": self.project_id,
            "project_root": str(self.project_root),
            "project_name": self.project_metadata.get("project_name"),
            "organism": self.project_metadata.get("organism"),
            "domain": self.project_metadata.get("domain"),
            "configs": str(self.configs),
            "data": str(self.data),
            "raw": str(self.raw),
            "external": str(self.external),
            "interim": str(self.interim),
            "processed": str(self.processed),
            "processed_rna": str(self.processed_rna),
            "processed_atac": str(self.processed_atac),
            "processed_multiome": str(self.processed_multiome),
            "processed_integration": str(self.processed_integration),
            "processed_grn": str(self.processed_grn),
            "processed_gnn": str(self.processed_gnn),
            "reports": str(self.reports),
            "tables": str(self.tables),
            "tables_rna": str(self.tables_rna),
            "tables_atac": str(self.tables_atac),
            "tables_annotation": str(self.tables_annotation),
            "figures": str(self.figures),
            "figures_rna": str(self.figures_rna),
            "figures_atac": str(self.figures_atac),
            "figures_annotation": str(self.figures_annotation),
            "manifests": str(self.manifests),
            "logs": str(self.logs),
            "models": str(self.models),
            "legacy_root": str(self.legacy_root),
            "legacy_data": str(self.legacy_data),
            "legacy_results": str(self.legacy_results),
        }  # χρήσιμο για manifests και debugging