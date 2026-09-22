"""Run test scenarios and load references; assertions live in pytest files."""

from pathlib import Path
import tomllib

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASELINES = ROOT / "tests" / "baselines"
CHECKS_PATH = ROOT / "tests" / "checks.toml"
with CHECKS_PATH.open("rb") as _file:
    CHECKS = tomllib.load(_file)


def run_scenario(name):
    """Run one configuration without a browser or plot interaction."""
    from examples.weathering import app

    overrides = CHECKS["scenarios"][name].get("overrides", {})
    _, definitions = app.run(defs=dict(overrides))
    return definitions


def load_baseline(name):
    path = BASELINES / f"{name}.npz"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing reference {path}. Tests never create references. "
            f"Review the model output, then run make update-baseline CASE={name}."
        )
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}
