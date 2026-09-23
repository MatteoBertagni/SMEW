PYTHON ?= .venv/bin/python
PYTEST_ARGS ?=
CASE ?=

.PHONY: build test test-serial example check-notebooks update-baseline

build:
	$(PYTHON) -m build

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
