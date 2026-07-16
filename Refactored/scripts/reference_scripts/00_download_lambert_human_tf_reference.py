# ============================================================
# 00_download_lambert_human_tf_reference.py
#
# Purpose
# -------
# Download and freeze the Lambert human transcription-factor
# reference files for reuse by downstream GRN scripts.
#
# Output directory
# ----------------
# Refactored/data/Human_TF_list/
#
# This script:
# - downloads the official Lambert Human TF v1.01 files,
# - validates that files are non-empty,
# - calculates SHA-256 checksums,
# - writes a YAML manifest,
# - writes a TSV checksum report.
#
# It does NOT:
# - load any PBMC data,
# - perform TF-to-RNA overlap,
# - apply expression filtering,
# - construct GRNs.
# ============================================================


# ============================================================
# 1. Imports
# ============================================================

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd
import yaml


# ============================================================
# 2. Resolve repository
# ============================================================

script_path = Path(__file__).resolve()

# Expected location:
#
# Refactored/
# └── scripts/
#     └── reference_scripts/
#         └── 00_download_lambert_human_tf_reference.py
#
# parents[0] -> reference_scripts
# parents[1] -> scripts
# parents[2] -> Refactored

refactor_dir = script_path.parents[2]

print("\nScript path:")
print(script_path)

print("\nRefactor repository:")
print(refactor_dir)

if not refactor_dir.is_dir():
    raise FileNotFoundError(
        "Refactor repository was not found:\n"
        f"{refactor_dir}"
    )


# ============================================================
# 3. Define output directory
# ============================================================

reference_dir = (
    refactor_dir
    / "data"
    / "Human_TF_list"
)

reference_dir.mkdir(
    parents=True,
    exist_ok=True,
)

print("\nReference output directory:")
print(reference_dir)


# ============================================================
# 4. Define official Lambert reference files
# ============================================================

reference_name = "Lambert human transcription factors"
reference_version = "v1.01"

tf_names_url = (
    "https://humantfs.ccbr.utoronto.ca/"
    "download/v_1.01/TF_names_v_1.01.txt"
)

tf_database_url = (
    "https://humantfs.ccbr.utoronto.ca/"
    "download/v_1.01/DatabaseExtract_v_1.01.csv"
)

tf_names_path = (
    reference_dir
    / "TF_names_v_1.01.txt"
)

tf_database_path = (
    reference_dir
    / "DatabaseExtract_v_1.01.csv"
)

manifest_path = (
    reference_dir
    / "lambert_human_tf_v1.01_manifest.yaml"
)

checksum_table_path = (
    reference_dir
    / "lambert_human_tf_v1.01_checksums.tsv"
)

run_metadata_path = (
    reference_dir
    / "lambert_human_tf_v1.01_download_metadata.json"
)


# ============================================================
# 5. Runtime configuration
# ============================================================

# False means that existing, non-empty files are retained.
overwrite_existing = False

download_timeout_seconds = 120

user_agent = (
    "MSc-PBMC-GRN-reference-downloader/1.0"
)


# ============================================================
# 6. Helper functions
# ============================================================

def utc_timestamp() -> str:
    """Return the current UTC time in ISO-8601 format."""
    return datetime.now(timezone.utc).isoformat()


def validate_file(
    path: Path,
    label: str,
) -> None:
    """Validate that a required file exists and is non-empty."""
    if not path.exists():
        raise FileNotFoundError(
            f"{label} was not found:\n"
            f"{path}"
        )

    if not path.is_file():
        raise ValueError(
            f"{label} is not a regular file:\n"
            f"{path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"{label} is empty:\n"
            f"{path}"
        )


def sha256_file(
    path: Path,
    chunk_size: int = 1024 * 1024,
) -> str:
    """Calculate SHA-256 without loading the whole file."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def download_file(
    url: str,
    output_path: Path,
    overwrite: bool = False,
) -> str:
    """
    Download a file safely using a temporary .part file.

    Returns
    -------
    str
        "downloaded" or "reused_existing"
    """
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if output_path.exists() and not overwrite:
        validate_file(
            output_path,
            label="Existing reference file",
        )

        print("\nExisting valid file retained:")
        print(output_path)
        print(
            f"Size: {output_path.stat().st_size:,} bytes"
        )

        return "reused_existing"

    temporary_path = output_path.with_suffix(
        output_path.suffix + ".part"
    )

    if temporary_path.exists():
        temporary_path.unlink()

    request = Request(
        url,
        headers={
            "User-Agent": user_agent,
        },
    )

    print("\nDownloading:")
    print(url)

    print("Destination:")
    print(output_path)

    try:
        with urlopen(
            request,
            timeout=download_timeout_seconds,
        ) as response:
            content = response.read()

        if not content:
            raise ValueError(
                "Remote server returned an empty response."
            )

        temporary_path.write_bytes(content)

        validate_file(
            temporary_path,
            label="Temporary downloaded file",
        )

        if output_path.exists():
            output_path.unlink()

        temporary_path.replace(
            output_path
        )

    except HTTPError as exc:
        if temporary_path.exists():
            temporary_path.unlink()

        raise RuntimeError(
            "HTTP error while downloading reference.\n"
            f"URL: {url}\n"
            f"HTTP status: {exc.code}\n"
            f"Reason: {exc.reason}"
        ) from exc

    except URLError as exc:
        if temporary_path.exists():
            temporary_path.unlink()

        raise RuntimeError(
            "Network error while downloading reference.\n"
            f"URL: {url}\n"
            f"Reason: {exc.reason}"
        ) from exc

    except Exception:
        if temporary_path.exists():
            temporary_path.unlink()

        raise

    validate_file(
        output_path,
        label="Downloaded reference file",
    )

    print("Download completed.")
    print(
        f"Size: {output_path.stat().st_size:,} bytes"
    )

    return "downloaded"


# ============================================================
# 7. Download official files
# ============================================================

started_at = utc_timestamp()

tf_names_status = download_file(
    url=tf_names_url,
    output_path=tf_names_path,
    overwrite=overwrite_existing,
)

tf_database_status = download_file(
    url=tf_database_url,
    output_path=tf_database_path,
    overwrite=overwrite_existing,
)


# ============================================================
# 8. Validate downloaded files
# ============================================================

validate_file(
    tf_names_path,
    label="Lambert TF names file",
)

validate_file(
    tf_database_path,
    label="Lambert complete database file",
)


# ============================================================
# 9. Basic content validation
# ============================================================

tf_names = pd.read_csv(
    tf_names_path,
    header=None,
    names=["tf_symbol"],
    dtype="string",
)

tf_names["tf_symbol"] = (
    tf_names["tf_symbol"]
    .str.strip()
)

if tf_names.empty:
    raise ValueError(
        "The Lambert TF-name table is empty."
    )

if tf_names["tf_symbol"].isna().any():
    raise ValueError(
        "The Lambert TF-name table contains missing symbols."
    )

if tf_names["tf_symbol"].eq("").any():
    raise ValueError(
        "The Lambert TF-name table contains empty symbols."
    )

duplicated_tf_symbols = (
    tf_names.loc[
        tf_names["tf_symbol"].duplicated(
            keep=False
        ),
        "tf_symbol",
    ]
    .dropna()
    .unique()
    .tolist()
)

if duplicated_tf_symbols:
    raise ValueError(
        "The Lambert TF-name table contains duplicated "
        "symbols.\n"
        f"Examples: {duplicated_tf_symbols[:20]}"
    )

tf_database = pd.read_csv(
    tf_database_path,
    low_memory=False,
)

if tf_database.empty:
    raise ValueError(
        "The Lambert complete database table is empty."
    )


# ============================================================
# 10. Calculate checksums
# ============================================================

tf_names_sha256 = sha256_file(
    tf_names_path
)

tf_database_sha256 = sha256_file(
    tf_database_path
)

print("\nReference checksums:")
print(f"- TF names: {tf_names_sha256}")
print(f"- Database: {tf_database_sha256}")


# ============================================================
# 11. Save checksum table
# ============================================================

checksum_table = pd.DataFrame(
    [
        {
            "reference_name": reference_name,
            "reference_version": reference_version,
            "file_role": "authoritative_tf_name_list",
            "filename": tf_names_path.name,
            "path": str(tf_names_path),
            "url": tf_names_url,
            "file_size_bytes": int(
                tf_names_path.stat().st_size
            ),
            "sha256": tf_names_sha256,
            "download_status": tf_names_status,
        },
        {
            "reference_name": reference_name,
            "reference_version": reference_version,
            "file_role": "complete_tf_database",
            "filename": tf_database_path.name,
            "path": str(tf_database_path),
            "url": tf_database_url,
            "file_size_bytes": int(
                tf_database_path.stat().st_size
            ),
            "sha256": tf_database_sha256,
            "download_status": tf_database_status,
        },
    ]
)

checksum_table.to_csv(
    checksum_table_path,
    sep="\t",
    index=False,
)


# ============================================================
# 12. Save YAML manifest
# ============================================================

manifest = {
    "reference_name": reference_name,
    "reference_version": reference_version,
    "reference_directory": str(reference_dir),
    "downloaded_or_validated_at_utc": utc_timestamp(),
    "authoritative_tf_catalogue": {
        "filename": tf_names_path.name,
        "path": str(tf_names_path),
        "url": tf_names_url,
        "sha256": tf_names_sha256,
        "file_size_bytes": int(
            tf_names_path.stat().st_size
        ),
        "n_tf_symbols": int(
            len(tf_names)
        ),
        "download_status": tf_names_status,
    },
    "complete_database": {
        "filename": tf_database_path.name,
        "path": str(tf_database_path),
        "url": tf_database_url,
        "sha256": tf_database_sha256,
        "file_size_bytes": int(
            tf_database_path.stat().st_size
        ),
        "shape": [
            int(tf_database.shape[0]),
            int(tf_database.shape[1]),
        ],
        "columns": [
            str(column)
            for column in tf_database.columns
        ],
        "download_status": tf_database_status,
    },
    "usage": {
        "authoritative_list_for_tf_overlap": (
            tf_names_path.name
        ),
        "notes": (
            "The complete database is retained for provenance "
            "and future TF annotation. Downstream PBMC scripts "
            "must not modify these frozen reference files."
        ),
    },
}

with manifest_path.open(
    "w",
    encoding="utf-8",
) as handle:
    yaml.safe_dump(
        manifest,
        handle,
        sort_keys=False,
        allow_unicode=True,
    )


# ============================================================
# 13. Save run metadata
# ============================================================

run_metadata = {
    "script": script_path.name,
    "script_path": str(script_path),
    "started_at_utc": started_at,
    "completed_at_utc": utc_timestamp(),
    "reference_directory": str(reference_dir),
    "overwrite_existing": overwrite_existing,
    "tf_names_status": tf_names_status,
    "tf_database_status": tf_database_status,
    "n_tf_symbols": int(len(tf_names)),
    "tf_database_shape": [
        int(tf_database.shape[0]),
        int(tf_database.shape[1]),
    ],
    "status": "completed",
}

with run_metadata_path.open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        run_metadata,
        handle,
        indent=2,
        ensure_ascii=False,
    )


# ============================================================
# 14. Validate all outputs
# ============================================================

required_outputs = {
    "Lambert TF names": tf_names_path,
    "Lambert complete database": tf_database_path,
    "Reference manifest": manifest_path,
    "Checksum table": checksum_table_path,
    "Run metadata": run_metadata_path,
}

print("\nValidating outputs:")

for label, path in required_outputs.items():
    validate_file(
        path,
        label=label,
    )

    print(f"- {label}: OK")
    print(f"  {path}")


# ============================================================
# 15. Final report
# ============================================================

print("\n" + "=" * 70)
print("LAMBERT HUMAN TF REFERENCE READY")
print("=" * 70)

print("\nReference directory:")
print(reference_dir)

print("\nReference version:")
print(reference_version)

print("\nCurated TF symbols:")
print(f"{len(tf_names):,}")

print("\nComplete database shape:")
print(
    f"{tf_database.shape[0]:,} rows × "
    f"{tf_database.shape[1]:,} columns"
)

print("\nFiles:")
print(f"- {tf_names_path.name}")
print(f"- {tf_database_path.name}")
print(f"- {manifest_path.name}")
print(f"- {checksum_table_path.name}")
print(f"- {run_metadata_path.name}")

print("\nNext step:")
print(
    "05A.1 will load TF_names_v_1.01.txt from this fixed "
    "reference directory and map it to the full-gene RNA "
    "object."
)