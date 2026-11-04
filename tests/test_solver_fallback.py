"""Exercise automatic fallback in fresh interpreters, before SMEW is imported."""

import os
import subprocess
import sys
import textwrap

import pytest


@pytest.mark.parametrize(("error_type", "module"), (
    ("ModuleNotFoundError", "smew._native._minpack"),
    ("ImportError", "smew._native._numba"),
    ("OSError", "smew._native._numba"),
    ("AttributeError", "smew._native._numba"),
    ("ValueError", "smew._native._numba"),
))
def test_native_load_failure_uses_scipy(error_type, module):
    script = textwrap.dedent("""
        import builtins
        import importlib.abc
        import logging
        import sys
        import warnings

        logging.basicConfig()
        import numpy as np
        from numba import config

        error_type, module = sys.argv[1:]
        class BrokenNative(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == module:
                    raise getattr(builtins, error_type)("simulated native load failure")
        sys.meta_path.insert(0, BrokenNative())

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            import smew
            from smew import _utils
            from smew.ic import f_CEC_Alk_to_conc
            from tests.test_ic import balance_inputs

            # Record the actual SciPy solves, including the nested CEC correction.
            solved = set()
            fsolve = _utils.fsolve
            def tracked_fsolve(equations, *args, **kwargs):
                solved.add(equations.__name__)
                return fsolve(equations, *args, **kwargs)
            _utils.fsolve = tracked_fsolve

            inputs = balance_inputs.__wrapped__()
            fractions, constants = smew.conc_to_f_CEC(inputs["conc_in"], 6., "loam", 1., 1.)
            np.testing.assert_allclose(fractions, inputs["f_CEC_in"], rtol=1e-9)
            result = smew.biogeochem_balance(**inputs)
            assert np.isfinite(result["pH"]).all()

            # Invert inventories to cover the 13-variable initialization solve.
            totals = [inputs["conc_in"][i] * (.4*.3*.6*1000) + fractions[i]*.1/charge
                      for i, charge in enumerate((2, 2, 1, 1))]
            _, initialized, _ = smew.total_to_f_CEC_and_conc(
                totals, 6., .07, [.6], "loam", .4, .3, .1, 1., 1.)
            np.testing.assert_allclose(initialized, fractions, rtol=1e-8)

            # The alkalinity correction uses its own scalar solver caller.
            conc = inputs["conc_in"]
            alkalinity = 2*conc[0] + 2*conc[1] + conc[2] + conc[3]
            corrected, _ = f_CEC_Alk_to_conc(fractions, constants, 6., alkalinity, 1., 1.)
            np.testing.assert_allclose(corrected, conc, rtol=1e-8)

        fallback = [w for w in caught if "Falling back to Python/SciPy fsolve" in str(w.message)]
        assert not fallback, caught
        assert solved == {"cec_calcium_equation", "water_alk_equations", "biogeochem_equations",
                          "total_to_cec_equations", "f_cec_alk_to_conc_equation"}, solved
        assert not config.DISABLE_JIT
        assert not hasattr(smew.biogeochem_balance_numba, "py_func")
        assert smew.f_CEC_to_conc.nopython_signatures
    """)
    process = subprocess.run(
        [sys.executable, "-c", script, error_type, module],
        env={**os.environ, "NUMBA_DISABLE_JIT": "0"},
        capture_output=True, text=True,
    )
    assert process.returncode == 0, process.stdout + process.stderr
    assert process.stderr.count("ERROR:smew._utils:SMEW native solver could not be loaded") == 1
    assert f"{error_type}: simulated native load failure" in process.stderr
    assert "Falling back to Python/SciPy fsolve" in process.stderr


def test_unrelated_native_error_is_not_hidden():
    script = textwrap.dedent("""
        import importlib.abc
        import sys

        class BrokenNative(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == "smew._native._numba":
                    raise RuntimeError("unrelated programming error")
        sys.meta_path.insert(0, BrokenNative())
        import smew
    """)
    process = subprocess.run(
        [sys.executable, "-c", script],
        env={**os.environ, "NUMBA_DISABLE_JIT": "0"},
        capture_output=True, text=True,
    )
    assert process.returncode != 0
    assert "RuntimeError: unrelated programming error" in process.stderr
    assert "Falling back" not in process.stderr
