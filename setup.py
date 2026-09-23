"""Build the private Cython equations extension from the maintained Python source."""

import os
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
    for suffix in (".py", ".pxd"):
        copyfile(ROOT / "smew" / f"equations{suffix}", generated / f"_equations{suffix}")

    extension = Extension(
        "smew._native._equations",
        [str(generated_relative / "_equations.py")],
    )
    return cythonize(
        [extension],
        include_path=[str(ROOT / "build" / "cython_sources")],
        compiler_directives={"language_level": 3, "boundscheck": False, "wraparound": False},
    )


extensions = native_extensions()

setup(
    ext_modules=extensions,
    options={"build": {"build_base": "build/smew_native" if extensions else "build/smew_python"}},
)
