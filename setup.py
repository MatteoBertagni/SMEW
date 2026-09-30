"""Build private Cython equations and the bundled cminpack solver."""

import os
import ast
import platform
from pathlib import Path
from shutil import copyfile

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.errors import OptionError


ROOT = Path(__file__).resolve().parent
CPU_TARGET = os.environ.get("SMEW_CPU_TARGET", "portable")
if CPU_TARGET not in {"portable", "native", "avx2", "avx512"}:
    raise ValueError("SMEW_CPU_TARGET must be portable, native, avx2, or avx512")


class BuildExt(build_ext):
    """Tune only explicitly requested local builds; published wheels stay portable."""

    def build_extensions(self):
        flags = []
        machine = platform.machine().lower()
        if CPU_TARGET != "portable":
            if CPU_TARGET in {"avx2", "avx512"} and machine not in {
                "x86_64", "amd64", "i386", "i686", "x86",
            }:
                raise OptionError(f"SMEW_CPU_TARGET={CPU_TARGET} requires an x86 target")
            if self.compiler.compiler_type == "msvc":
                if CPU_TARGET == "native":
                    raise OptionError(
                        "MSVC has no native CPU auto-detection flag. Set "
                        "SMEW_CPU_TARGET=avx2 or avx512 only if your CPU supports "
                        "it, or use portable."
                    )
                flags = [f"/arch:{CPU_TARGET.upper()}"]
            elif self.compiler.compiler_type in {"unix", "mingw32"}:
                if CPU_TARGET == "native":
                    flags = ["-mcpu=native" if machine in {"arm64", "aarch64"}
                             else "-march=native"]
                else:
                    flags = ["-mavx2" if CPU_TARGET == "avx2" else "-mavx512f"]
            else:
                raise OptionError(
                    f"CPU tuning is unsupported with {self.compiler.compiler_type}"
                )
        for extension in self.extensions:
            extension.extra_compile_args = [*extension.extra_compile_args, *flags]
        # Distutils does not track compiler flags when deciding to reuse objects.
        # In particular, never reuse tuned objects in a later portable build.
        self.force = True
        super().build_extensions()


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
    source = (ROOT / "smew" / "equations.py").read_text(encoding="utf-8")
    # SciPy's vector wrappers pass NumPy arrays; the native residuals instead
    # receive C pointers. Keep the scientific functions verbatim in both builds.
    python_only = {"biogeochem_equations", "total_to_cec_equations", "kelland_equations"}
    lines = source.splitlines(keepends=True)
    for node in ast.parse(source).body:
        if isinstance(node, ast.FunctionDef) and node.name in python_only:
            for line in range(node.lineno - 1, node.end_lineno):
                lines[line] = "\n"
    (generated / "_equations.py").write_text("".join(lines), encoding="utf-8")
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
        # Numba resolves this symbol with ctypes.CDLL, including on Windows.
        export_symbols=["smew_solve"],
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
    cmdclass={"build_ext": BuildExt},
    options={"build": {"build_base": "build/smew_native" if extensions else "build/smew_python"}},
)
