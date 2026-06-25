from pathlib import Path
import json
import pandas as pd


def read_csv(path, **kwargs):
    """
    Read CSV with path validation.
    """
    
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"CSV file not found: {path}")
    
    return pd.read_csv(path, **kwargs)


def write_csv(df, path, index=False, **kwargs):
    """
    Write DataFrame to CSV and create parent directories.
    """
    
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    df.to_csv(path, index=index, **kwargs)


def read_json(path):
    """
    Read JSON file.
    """
    
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"JSON file not found: {path}")
    
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(obj, path, indent=2):
    """
    Write JSON file and create parent directories.
    """
    
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=indent, ensure_ascii=False)


def write_manifest(manifest, path):
    """
    Write a run/output manifest as JSON.
    """
    
    write_json(manifest, path, indent=2)
