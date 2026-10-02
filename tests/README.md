# Model tests

From the repository root, using Python 3.11 or newer:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
make test
```

## Installed wheel checks

The packaging workflow installs each wheel in a fresh environment and copies
the tests to a temporary directory outside the checkout. It runs:

```bash
python tests/check_wheel.py
python -m pytest --import-mode=importlib tests/test_native_equations.py tests/test_native_solver.py tests/test_backend.py
```

`check_wheel.py` requires both extensions, their Numba binding, and the license
notices to be installed from the wheel. It rejects editable installations. The
pytest checks exercise native equations, solver buffers, compiled initialization,
and short simulations against the Python backend. This subset needs pytest and
marimo, but not the full-year baseline files. Source archives include this test
subset's supporting files; full regressions require a repository checkout.

## Full model regressions

Pytest runs three configurations of the [marimo example](../tests/example_marimo_notebook.py):
the defaults, constant moisture, and no added rock. Each covers **365 days at a
ten-minute timestep** (52,560 samples). Fixed weekly rainfall makes the runs
repeatable. `make test` runs the suite in separate processes with `NUMBA_DISABLE_JIT=1`
and `NUMBA_DISABLE_JIT=0`, using the same references and tolerances. Compiled regression
tests fail if the native extensions are unavailable; they are not skipped.
Assertions compare every selected timestep with its reference and
check finite values and physical ranges separately. Solver warnings stay visible.
`make test` runs Python mode in parallel, selecting the number of workers
automatically and capping it at the three scenarios. Tests are grouped by
scenario, ensuring each expensive model configuration executes once per backend.
Compiled mode runs in one process so the scenarios and backend probe share
Numba compilations. Its process starts with `NUMBA_NRT_STATS=1` for allocation
checks. Use `make test-serial` to also run Python mode without multiprocessing.

```bash
make test PYTEST_ARGS='-k physical'  # range checks without references
make test-serial                    # disable multiprocessing
make test PYTHON=python             # use an activated environment
make check-notebooks                # validate the marimo notebook
make example                       # open the example; click to generate plots
```

The example explains its inputs directly in Python. It can also run as a script
or be imported: `_, definitions = app.run()` exposes `definitions["results"]`.
Plotting stays inactive in tests.

## Configuration

[checks.toml](checks.toml) is the only TOML configuration file. It contains the
three scenarios, selected variables, units, tolerances and bounds. Comparisons use
`abs(current-reference) <= atol + rtol*abs(reference)`; pH uses absolute tolerance.
`bound_atol` separately allows small roundoff around physical bounds. Alkalinity
may be negative. Temperature and pH ranges are envelopes for these examples.
New scenarios belong under `[scenarios]`; new selected variables need a
`[variables]` rule and an entry in the notebook's `results`.

The numerical probes in `test_backend.py` compare helpers, initialization, and
short simulations with `rtol=1e-9` and `atol=1e-9` (float32 seasonal vegetation
uses `rtol=1e-6`). Python mode runs in a fresh process; compiled mode reuses the
pytest process when JIT and NRT statistics are enabled at startup. Other direct
pytest invocations use a fresh compiled process to retain the allocation checks.
A separate fresh process checks that compiled helpers need no native solver.
The chemistry probes check that empty and nonempty pore arrays share the same
compiled calculation. The Python entry point converts missing pores to empty
arrays; the public compiled entry point requires arrays directly.
The probes run once from the compiled test suite and exercise both startup
modes. This allows small accumulated numerical differences between
SciPy/MINPACK and Numba/cminpack across platforms; a relative tolerance of
`1e-10` proved too tight in CI.

The compiled probe uses `NUMBA_NRT_STATS` to check that repeated
successful and failed calls leave the number of live Numba allocations unchanged.
It covers native residual failures, insufficient cations, missing nonlinear
pore inputs, zero-area nonlinear distributions, and inactive-soil respiration errors. These
checks run after compilation warmup; compiler memory and process RSS are not
used as a proxy for leaked arrays.

`test_errors.py` checks code uniqueness, generated messages, exception
serialization, and propagation through a parent Numba
function. Backend probes also assert the specific codes produced by failed
simulations while checking that their arrays are released.

## Updating references

`make test` never writes references (baselines); missing references fail.
For an intentional change, review the reported differences, then explicitly
update each affected reference:

```bash
make update-baseline CASE=default
make test
```

This checks physical ranges, reports maximum changes, and replaces the selected
`baselines/<scenario>.npz`. These `.npz` files contain arrays only; there are no
per-scenario TOML metadata files. Review and commit changed references with the
reason for the change; Git retains their previous versions.
