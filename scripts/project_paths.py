"""Paths for local, non-public inputs and generated artifacts.

No proprietary data are distributed with this repository.
Set FANKOU_DATA_DIR and optionally FANKOU_OUTPUT_DIR to use your own files.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("FANKOU_DATA_DIR", str(Path.home() / "fankou_private_inputs"))).expanduser()
OUTPUT_DIR = Path(os.environ.get("FANKOU_OUTPUT_DIR", str(REPO_ROOT / "outputs"))).expanduser()

def required_file(name: str) -> str:
    path = DATA_DIR / name
    if not path.is_file():
        raise FileNotFoundError(
            f"Required private input is missing: {path}. "
            "Set FANKOU_DATA_DIR to the folder containing your confidential inputs (see README.md)."
        )
    return str(path)

def output_path(name: str) -> str:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return str(OUTPUT_DIR / name)
