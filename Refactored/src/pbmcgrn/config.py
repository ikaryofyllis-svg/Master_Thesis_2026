from pathlib import Path
import yaml


def load_yaml(path):
    """
    Load a YAML config file.
    """
    
    path = Path(path)
    
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def save_yaml(config, path):
    """
    Save a dictionary as YAML.
    """
    
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(
            config,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )


def load_project_config(refactor_root):
    """
    Load configs/project.yaml from the refactored project.
    """
    
    refactor_root = Path(refactor_root)
    config_path = refactor_root / "configs" / "project.yaml"
    
    return load_yaml(config_path)


def load_cell_type_config(refactor_root):
    """
    Load configs/cell_types.yaml from the refactored project.
    """
    
    refactor_root = Path(refactor_root)
    config_path = refactor_root / "configs" / "cell_types.yaml"
    
    return load_yaml(config_path)
