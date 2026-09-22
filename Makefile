PYTHON ?= .venv/bin/python
PYTEST_ARGS ?=
CASE ?=

.PHONY: test test-serial example check-notebooks update-baseline

test:
	$(PYTHON) -m pytest -n auto --maxprocesses=3 --dist=loadgroup tests $(PYTEST_ARGS)

test-serial:
	$(PYTHON) -m pytest tests $(PYTEST_ARGS)

example:
	$(PYTHON) -m marimo edit examples/weathering.py

check-notebooks:
	$(PYTHON) -m marimo check --strict examples/weathering.py

# Explicitly select one configured scenario; `make test` never updates references.
update-baseline:
	$(PYTHON) -m tests.update_baseline --case "$(CASE)"
