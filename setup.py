from setuptools import setup, Extension
import os
import sys
import logging
import warnings
from setuptools import setup, Extension

# Set up a logger that works within the setuptools ecosystem
logger = logging.getLogger(__name__)

# If we are running as a script, we need a basic config, 
# otherwise setuptools/pip handles the handlers.
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def check_header(path):
    """Helper to verify the header actually exists at the given path."""
    exists = path and os.path.exists(os.path.join(path, "cminpack.h"))
    if exists:
        logger.debug(f"Confirmed cminpack.h at: {path}")
    return exists

# Allow overriding paths via environment variables
cminpack_include = os.environ.get("CMINPACK_INCLUDE")
cminpack_lib = os.environ.get("CMINPACK_LIB")

# Automatically detect Conda/Virtualenv paths
if not check_header(cminpack_include):
    search_paths = [
        os.path.join(sys.prefix, "include", "cminpack-1"),
        os.path.join(sys.prefix, "include"),
        "/opt/homebrew/include/cminpack-1",  # Apple Silicon Homebrew
        "/usr/local/include/cminpack-1",     # Intel Homebrew
        "/usr/include/cminpack-1",           # Ubuntu/Debian
        "/usr/include",                      # Generic Linux fallback
        "/usr/local/include",
    ]

    cminpack_include = None
    for path in search_paths:
        if check_header(path):
            cminpack_include = path
            break

# If headers are found, determine lib path conditionally
if cminpack_include and not cminpack_lib:
    search_paths = [
        os.path.join(sys.prefix, "lib"),
        "/opt/homebrew/lib",
        "/usr/local/lib",
        "/usr/lib",
    ]
    for path in search_paths:
        # Check for standard shared/static library file extensions
        if any(os.path.exists(os.path.join(path, f)) for f in ["libcminpack.so", "libcminpack.dylib", "cminpack.lib", "libcminpack.a"]):
            cminpack_lib = path
            break

extensions = []

# Only attempt to build the Cython extension if cminpack headers were definitively found
if cminpack_include:
    logger.info(f"cminpack found. Headers: {cminpack_include}, Libs: {cminpack_lib}")
    
    try:
        import numpy
        from Cython.Build import cythonize
    except ImportError as e:
        missing_pkg = "numpy" if "numpy" in str(e) else "Cython"
        logger.warn(f"{missing_pkg} not installed. Skipping fast solver build.")
    else:
        try:
            include_dirs = [numpy.get_include(), cminpack_include]
            library_dirs = [cminpack_lib] if cminpack_lib else []
            runtime_library_dirs = [cminpack_lib] if cminpack_lib else []

            ext = Extension(
                name="smew.biogeochem_fast",
                sources=["smew/biogeochem_fast.pyx"],
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                libraries=["cminpack"],
                runtime_library_dirs=runtime_library_dirs,
                extra_compile_args=['-O3', '-mavx2']
            )
        
            extensions = cythonize([ext], language_level="3")
        except Exception as e:
            logger.error(f"Failed to cythonize extension: {e}")
            extensions = []
else:
    logger.warning("cminpack headers not found. Building SMEW without the fast Cython solver.")

setup(
    ext_modules=extensions
)