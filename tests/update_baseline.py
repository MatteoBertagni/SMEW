"""Explicitly replace the numerical reference for one configured scenario."""

import argparse

import numpy as np

from tests.support import BASELINES, CHECKS, ROOT, run_scenario
from tests.test_physical_ranges import assert_physical_range


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=tuple(CHECKS["scenarios"]))
    args = parser.parse_args()

    print(f"Running {args.case} and checking physical ranges.", flush=True)
    results = run_scenario(args.case)["results"]
    arrays = {"time_days": results["time_days"]}

    for variable, rule in CHECKS["variables"].items():
        values = results[variable]
        assert_physical_range(values, rule, label=f"{args.case}/{variable}")
        if values.ndim < 1 or values.shape[-1] != len(arrays["time_days"]):
            raise ValueError(f"{variable}: last axis must correspond to time")
        arrays[variable] = values

    path = BASELINES / f"{args.case}.npz"
    if path.exists():
        print("Maximum absolute changes from the previous reference:")
        with np.load(path, allow_pickle=False) as old:
            for variable, values in arrays.items():
                if variable in old and old[variable].shape == values.shape:
                    change = np.max(np.abs(values - old[variable]))
                    print(f"  {variable}: {change:.8g}")
                else:
                    print(f"  {variable}: new output or changed shape")

    BASELINES.mkdir(exist_ok=True)
    np.savez_compressed(path, **arrays)
    print(f"Wrote {path.relative_to(ROOT)}.")
    print("Review the change and run make test before committing it.")


if __name__ == "__main__":
    main()
