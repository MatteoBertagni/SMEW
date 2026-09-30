"""Fail wheel CI on missing binaries/licenses or imports from a source checkout.

Run from outside the checkout, before pytest (whose native tests may skip when
extensions are absent). The numerical and Numba integration checks are reused
from test_native_equations.py, test_native_solver.py, and test_backend.py.
"""

from importlib.metadata import distribution
from pathlib import Path
import json

import smew
import pyeto
from smew._native import _equations, _minpack
from smew._simulation_compiled import _native_solve


installed = distribution("smew")
direct_url = json.loads(installed.read_text("direct_url.json") or "{}")
assert not direct_url.get("dir_info", {}).get("editable"), "Testing an editable install"
files = installed.files
assert files, "Missing wheel RECORD"
recorded = {Path(installed.locate_file(item)).resolve() for item in files}
for module in (smew, pyeto, _equations, _minpack):
    assert Path(module.__file__).resolve() in recorded, module.__file__
assert _native_solve is not None
assert any(str(item).endswith("/licenses/LICENSE.txt") for item in files)
assert any(str(item).endswith("/licenses/third_party/cminpack/CopyrightMINPACK.txt")
           for item in files)
assert not any(str(item).endswith((".nbc", ".nbi")) for item in files)
print(f"Testing installed SMEW {installed.version}: {smew.__file__}")
