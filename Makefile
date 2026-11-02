PYTHON ?= python
PYTEST_ARGS ?=
CASE ?=
SPHINX_ARGS ?=

.PHONY: build build-native build-python install-native install-python clean clean-build clean-docs test test-serial example check-notebooks update-baseline docs docs-check

# Build the standard SMEW packages for distribution.
build: build-native

# Build packages with the C extensions for faster model calculations.
build-native:
	SMEW_BUILD_NATIVE=1 $(PYTHON) -m build

# Build packages without C extensions for a Python-only installation.
build-python:
	SMEW_BUILD_NATIVE=0 $(PYTHON) -m build

# Rebuild and install C extensions after changing equations or solver code.
install-native:
	SMEW_BUILD_NATIVE=1 $(PYTHON) -m pip install --no-cache-dir -e .

# Install editable Python code without building C extensions.
install-python:
	SMEW_BUILD_NATIVE=0 $(PYTHON) -m pip install --no-cache-dir -e .

# Remove generated package and documentation files for a fresh start.
clean: clean-build clean-docs

# Remove generated packages and C extensions before a fresh build.
clean-build:
	rm -rf -- build dist smew.egg-info
	rm -f -- smew/_native/_equations*.so smew/_native/_equations*.pyd smew/_native/_minpack*.so smew/_native/_minpack*.pyd smew/_native/_minpack.c

# Remove generated documentation and cached notebook results before rebuilding.
clean-docs:
	rm -rf -- docs/_build docs/_generated docs/.jupyter_cache

# Check both model backends, using parallel workers to speed up Python tests.
test:
	NUMBA_DISABLE_JIT=1 $(PYTHON) -m pytest -n auto --maxprocesses=3 --dist=loadgroup tests --ignore=tests/test_backend.py $(PYTEST_ARGS)
	NUMBA_DISABLE_JIT=0 NUMBA_NRT_STATS=1 $(PYTHON) -m pytest -n 0 tests $(PYTEST_ARGS)

# Check both model backends without parallel workers to reduce memory use.
test-serial:
	NUMBA_DISABLE_JIT=1 $(PYTHON) -m pytest tests --ignore=tests/test_backend.py $(PYTEST_ARGS)
	NUMBA_DISABLE_JIT=0 NUMBA_NRT_STATS=1 $(PYTHON) -m pytest tests $(PYTEST_ARGS)

# Open the interactive example to explore model settings and results.
example:
	$(PYTHON) -m marimo edit tests/example_marimo_notebook.py

# Check the marimo example for structural errors without running calculations.
check-notebooks:
	$(PYTHON) -m marimo check --strict tests/example_marimo_notebook.py

# Build the documentation website to preview pages and notebook results.
docs:
	$(PYTHON) -m sphinx -b html -W --keep-going $(SPHINX_ARGS) docs docs/_build/html

# Rebuild all pages and rerun the example to catch documentation errors.
docs-check:
	$(PYTHON) -m sphinx -b html -W --keep-going -E -a -D nb_execution_mode=force $(SPHINX_ARGS) docs docs/_build/html

# Replace reference results for CASE after reviewing an intentional model change.
update-baseline:
	$(PYTHON) -m tests.update_baseline --case "$(CASE)"
