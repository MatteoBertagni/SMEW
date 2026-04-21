from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os

# Define the C-extension
extensions = [
    Extension(
        name="smew.biogeochem_fast",  # Places the compiled .so file inside the smew folder
        sources=["smew/biogeochem_fast.pyx"],
        include_dirs=[
            numpy.get_include(),
            "/usr/include/cminpack-1" # Default Ubuntu path. If using conda, you may need to adjust this.
        ],
        libraries=["cminpack"], # Tells the linker to look for libcminpack.so
    )
]

setup(
    ext_modules=cythonize(extensions, language_level="3")
)
