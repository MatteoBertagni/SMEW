# Model failures

TODO: build proper sphinx doc

`ErrorCode` values are globally unique and grouped by component: 1000–1999
for biogeochemistry, 2000–2999 for respiration, and 3000–3999 for initialization
solvers. Zero means success. Published codes must not be renumbered or reused.
Compare enum members in code rather than parsing error messages.

The public error API is `smew.ErrorCode`, `smew.SMEWError`, and
`smew.raise_for_error(code, message)`. Each failure branch constructs its own
complete debug message inside the compiled calculation.

## Python and Dask

Python entry points validate essential inputs before calling compiled code.
Each component has one `*_INVALID_INPUT` code; its message describes the input
problem and relevant names or shapes. Numerical failures raise `SMEWError`
subclasses after compiled code returns normally. Existing `except ValueError`
and `except RuntimeError` handlers still catch the corresponding failures.

```python
import smew

try:
    result = smew.biogeochem_balance(**inputs)
except smew.SMEWError as exc:
    print(exc.code)       # ErrorCode member; int(exc.code) gives its number
    print(exc.code.name)  # e.g. CHEMISTRY_RESIDUAL
    print(exc.message)    # complete debug message returned by compiled code
    print(str(exc))       # same message prefixed by the code and its name
    raise
```

For example:

```text
[SMEW1121 CHEMISTRY_RESIDUAL] Residual evaluation failed (nonfinite values); solver_status=-1, max_abs_residual=nan; biogeochem.chemistry; step=1, time_days=0.041666667, previous_pH=4, trial_H=100, moisture=0.49987647, temperature_C=3.8210268
```

`previous_pH` describes the preceding completed timestep. `trial_H` is from the
failed solver attempt, not an accepted state. Timestep indices are zero-based;
aggregate and input failures omit timestep information. Native solver flags
appear in the message separately from the SMEW error code. Residual norms are
omitted when buffers may be unwritten and reported as `nan` for nonfinite
residuals.

The internal `_float_text` helper formats debug floats to approximately eight
significant digits, including scientific notation, NaN and infinity, entirely
in Numba. Each failure can concatenate as many values as needed.

Exceptions preserve their code and message through pickle for transport by
Dask. The error payload contains no simulation arrays. Python tracebacks may
still hold frame locals, so persistent sweep logs should store the code and
message rather than exception traceback objects.

## Calling from another Numba model

The public compiled entry points are:

- `smew.biogeochem_balance_numba`
- `smew.respiration_numba`
- `smew.conc_to_f_CEC_numba`
- `smew.total_to_f_CEC_and_conc_numba`
- `smew.Kelland_numba`

Every entry point returns `(result, error_code, error_message)`:

- Failure: `result=None`, a nonzero integer code, and a complete message.
- Success: the numerical result and code `0`. The message is normally empty.
  Initialization solvers may return a nonfatal MINPACK warning message with
  code `0`; Python wrappers emit it as `RuntimeWarning`. Compiled callers can
  propagate that message to their own Python boundary.

Use the code to distinguish failure from success. Successful biogeochemistry
results are `(name, value)` pairs, converted to a dictionary by its Python
wrapper. Other entry points preserve their numerical result tuples/lists.
Partial trajectories are not returned on failure.

The compiled entry points assume inputs have already been validated by their
caller. They report numerical failures during computation and do not repeat
the Python input checks. When integrating with another compiled model, validate
inputs before entering that model and preserve these requirements:

- Array arguments are one-dimensional contiguous `float64` arrays.
- Biogeochemistry time series are nonempty and have the same length. Initial
  concentrations and exchange constants contain 5 values each; exchange
  fractions contain 6. `dt` is finite and positive.
- With rock added, `0 <= t_app / dt < len(s)`. Mineral names are supported,
  supplied as a nonempty tuple, and match the rock fraction array's length.
  Particle diameters are nonempty and match the particle fraction array's
  length. For no rock, use `("",)` and empty arrays for unused particle inputs.
- The surface model is `constant`, `linear`, or `nonlinear`. Nonlinear scaling
  requires matching, nonempty pore arrays and `0 < mixalf_in <= 1`. Pore arrays
  may be `None` for constant/linear scaling.
- Respiration uses nonempty moisture, vegetation, and temperature arrays of
  equal length, with a supported soil name.
- Initialization solvers use 5 concentrations, 4 totals where applicable,
  nonempty moisture arrays where applicable, and a soil supported by the CEC
  constants (`smew.constants.CEC_SOIL_TYPES`).

Violations of these input requirements are outside the compiled error contract;
they may cause compilation errors, invalid results, or out-of-bounds access.

```python
from numba import njit
import smew

@njit
def parent_model(s, vegetation, soil_temperature):
    result, code, message = smew.respiration_numba(
        None, 1000.0, None, 1.0, "loam", s, vegetation,
        1.0, 0.3, soil_temperature, 1.0, 1.0, 1000.0,
    )
    if code != smew.ErrorCode.OK or result is None:
        return None, code, message
    soc, microbial, root, diffusivity = result
    return soc, code, message

soc, code, message = parent_model(s, vegetation, soil_temperature)
smew.raise_for_error(code, message)  # Python boundary only
```

Compiled parents should propagate failures by returning normally. Raising while
a compiled parent owns arrays can reintroduce Numba's exception allocation leak.

This contract applies to the entry points listed above. Standalone helpers and
the legacy two-application model retain their existing exception interfaces.
Unexpected programming errors, compilation errors, and process termination are
not converted into numerical failure returns.
