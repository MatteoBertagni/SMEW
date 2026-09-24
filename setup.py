"""Build private Cython equations and the bundled cminpack solver."""

import os
import ast
from pathlib import Path
from shutil import copyfile

from setuptools import Extension, setup


ROOT = Path(__file__).resolve().parent


def native_extensions():
    setting = os.environ.get("SMEW_BUILD_NATIVE", "1")
    if setting == "0":
        return []
    if setting != "1":
        raise ValueError("SMEW_BUILD_NATIVE must be '0' or '1'")

    from Cython.Build import cythonize

    generated_relative = Path("build") / "cython_sources" / "smew" / "_native"
    generated = ROOT / generated_relative
    generated.mkdir(parents=True, exist_ok=True)
    source = (ROOT / "smew" / "equations.py").read_text()
    # SciPy's vector wrappers pass NumPy arrays; the native residuals instead
    # receive C pointers. Keep the scientific functions verbatim in both builds.
    python_only = {"biogeochem_equations", "total_to_cec_equations", "kelland_equations"}
    lines = source.splitlines(keepends=True)
    for node in ast.parse(source).body:
        if isinstance(node, ast.FunctionDef) and node.name in python_only:
            for line in range(node.lineno - 1, node.end_lineno):
                lines[line] = "\n"
    (generated / "_equations.py").write_text("".join(lines))
    copyfile(ROOT / "smew" / "equations.pxd", generated / "_equations.pxd")

    equations = Extension(
        "smew._native._equations",
        [str(generated_relative / "_equations.py")],
    )
    cminpack = Path("third_party") / "cminpack"
    solver = Extension(
        "smew._native._minpack",
        ["smew/_native/_minpack.pyx"] + [
            str(cminpack / f"{name}.c") for name in (
                "hybrd", "dogleg", "dpmpar", "enorm", "fdjac1",
                "qform", "qrfac", "r1mpyq", "r1updt",
            )
        ],
        include_dirs=[str(cminpack)],
        define_macros=[
            ("CMINPACK_NO_DLL", "1"), ("__cminpack_double__", "1"),
        ],
        libraries=["m"] if os.name == "posix" else [],
    )
    return cythonize(
        [equations, solver],
        build_dir="build/cython_generated",
        include_path=[str(ROOT / "build" / "cython_sources")],
        compiler_directives={"language_level": 3, "boundscheck": False,
                             "wraparound": False, "cdivision": True,
                             "cpow": True},
    )


extensions = native_extensions()

setup(
    ext_modules=extensions,
    options={"build": {"build_base": "build/smew_native" if extensions else "build/smew_python"}},
)
