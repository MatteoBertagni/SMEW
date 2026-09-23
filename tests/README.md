# Model tests

From the repository root, using Python 3.11 or newer:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
make test
```

Pytest runs three configurations of the [marimo example](../tests/example_marimo_notebook.py):
the defaults, constant moisture, and no added rock. Each covers **365 days at a
ten-minute timestep** (52,560 samples). Fixed weekly rainfall makes the runs
repeatable. Assertions compare every selected timestep with its reference and
check finite values and physical ranges separately. Solver warnings stay visible.
`make test` selects the number of workers automatically and caps it at the three
scenarios, so a two-core runner uses two workers. Tests are grouped by scenario,
ensuring each expensive model configuration executes once. Use `make test-serial`
to run the tests one by one for debugging or memory-constrained machines.

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
