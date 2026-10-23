# SPDX-License-Identifier: AGPL-3.0-only
"""Model error codes and Python exceptions.

Compiled entry points return (result, error_code, error_message).
Never renumber or reuse a published code. Raise only at the Python boundary.
"""

from enum import IntEnum, unique

import numpy as np
from numba.extending import register_jitable


@unique
class ErrorCode(IntEnum):
    OK = 0
    # Each nonlinear solve has its own codes, including the fallback attempt.
    RAIN_RESIDUAL = 1101
    RAIN_SOLVER_INPUT = 1102
    HYDROGEN_RESIDUAL = 1111
    HYDROGEN_SOLVER_INPUT = 1112
    CHEMISTRY_RESIDUAL = 1121
    CHEMISTRY_SOLVER_INPUT = 1122
    CHEMISTRY_RETRY_RESIDUAL = 1131
    CHEMISTRY_RETRY_SOLVER_INPUT = 1132
    CHEMISTRY_NO_CONVERGENCE = 1133
    BIOGEOCHEM_INSUFFICIENT_CATIONS = 1201
    BIOGEOCHEM_MISSING_PORES = 1204
    # Organic carbon: 2000-2999.
    RESPIRATION_INITIAL_ACTIVITY = 2001
    RESPIRATION_MEAN_ACTIVITY = 2002
    # Initial condition solvers: 3000-3999.
    CEC_RESIDUAL = 3001
    CEC_SOLVER_INPUT = 3002
    TOTAL_CEC_RESIDUAL = 3011
    TOTAL_CEC_SOLVER_INPUT = 3012
    KELLAND_RESIDUAL = 3021
    KELLAND_SOLVER_INPUT = 3022
    WEATHERING_ERROR = 4001


@register_jitable
def _float_text(value):
    """Format debug floats to about eight significant digits in nopython mode.

    Numba cannot format floats with str/f-strings. Scaling a decade at a time
    also handles subnormal values without overflowing powers of ten.
    """
    if np.isnan(value):
        return "nan"
    negative = np.signbit(value)
    sign = "-" if negative else ""
    if np.isinf(value):
        return sign + "inf"
    if value == 0:
        return sign + "0"
    magnitude = abs(float(value))
    exponent = 0
    while magnitude >= 10.0:
        magnitude /= 10.0
        exponent += 1
    while magnitude < 1.0:
        magnitude *= 10.0
        exponent -= 1
    rounded = int(round(magnitude * 10000000.0))
    if rounded >= 100000000:
        rounded = 10000000
        exponent += 1
    digits = str(rounded).rstrip("0")
    if exponent < -4 or exponent >= 8:
        mantissa = digits[0]
        if len(digits) > 1:
            mantissa += "." + digits[1:]
        return sign + mantissa + "e" + ("+" if exponent >= 0 else "-") + str(abs(exponent))
    point = exponent + 1
    if point <= 0:
        return sign + "0." + "0" * (-point) + digits
    if point >= len(digits):
        return sign + digits + "0" * (point - len(digits))
    return sign + digits[:point] + "." + digits[point:]


@register_jitable
def _solver_error(status, residual, residual_code, input_code, context):
    """Return a failure triple without reading unwritten solver buffers."""
    code = residual_code if status == -1 else input_code
    reason = ("Residual evaluation failed (nonfinite values)" if status == -1
              else "Native solver received invalid buffers or settings")
    message = reason + "; solver_status=" + str(status)
    if status == -1 or status > 0:
        norm = 0.0
        for value in residual:
            if not np.isfinite(value):
                norm = np.nan
                break
            norm = max(norm, abs(value))
        message += ", max_abs_residual=" + _float_text(norm)
    return None, code.value, message + "; " + context


class SMEWError(Exception):
    """Model failure with a globally unique code and a complete debug message."""

    def __init__(self, code, message):
        self.code = ErrorCode(code)
        self.message = message
        super().__init__(f"[SMEW{int(self.code):04d} {self.code.name}] {message}")

    def __reduce__(self):
        # Preserve constructor arguments when Dask/pickle transports the error.
        return type(self), (self.code, self.message)


class SMEWValueError(SMEWError, ValueError):
    """Model input or convergence error; remains catchable as ValueError."""


class SMEWSolverError(SMEWError, RuntimeError):
    """Invalid solver residual; remains catchable as RuntimeError."""


def raise_for_error(code, message):
    """Raise a Python exception after the compiled calculation has returned."""
    if code != ErrorCode.OK:
        exception = SMEWSolverError if code in (
            ErrorCode.RAIN_RESIDUAL, ErrorCode.HYDROGEN_RESIDUAL,
            ErrorCode.CHEMISTRY_RESIDUAL, ErrorCode.CHEMISTRY_RETRY_RESIDUAL,
            ErrorCode.CEC_RESIDUAL, ErrorCode.TOTAL_CEC_RESIDUAL,
            ErrorCode.KELLAND_RESIDUAL,
        ) else SMEWValueError
        raise exception(code, message)
