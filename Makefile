PYTHON ?= .venv/bin/python
PYTEST_ARGS ?=
CASE ?=

.PHONY: build build-native build-python install-native install-python clean test test-serial example check-notebooks update-baseline

build: build-native

build-native:
	SMEW_BUILD_NATIVE=1 $(PYTHON) -m build

build-python:
	SMEW_BUILD_NATIVE=0 $(PYTHON) -m build

install-native:
	SMEW_BUILD_NATIVE=1 $(PYTHON) -m pip install --no-cache-dir -e .

install-python:
	SMEW_BUILD_NATIVE=0 $(PYTHON) -m pip install --no-cache-dir -e .

clean:
	rm -rf -- build dist smew.egg-info
	rm -f -- smew/_native/_equations*.so smew/_native/_equations*.pyd smew/_native/_minpack*.so smew/_native/_minpack*.pyd smew/_native/_minpack.c

test:
	$(PYTHON) -m pytest -n auto --maxprocesses=3 --dist=loadgroup tests $(PYTEST_ARGS)

test-serial:
	$(PYTHON) -m pytest tests $(PYTEST_ARGS)

example:
	$(PYTHON) -m marimo edit tests/example_marimo_notebook.py

check-notebooks:
	$(PYTHON) -m marimo check --strict tests/example_marimo_notebook.py

# Explicitly select one configured scenario; `make test` never updates references.
update-baseline:
	$(PYTHON) -m tests.update_baseline --case "$(CASE)"
