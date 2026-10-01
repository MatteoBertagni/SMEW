import pytest

from tests.support import CHECKS, load_baseline, run_scenario


@pytest.fixture(
    scope="session",
    params=[
        pytest.param(name, marks=pytest.mark.xdist_group(name))
        for name in CHECKS["scenarios"]
    ],
)
def scenario_name(request):
    return request.param


@pytest.fixture(scope="session")
def notebook_run(scenario_name):
    return run_scenario(scenario_name)


@pytest.fixture(scope="session")
def results(notebook_run):
    return notebook_run["results"]


@pytest.fixture(scope="session")
def baseline(scenario_name):
    return load_baseline(scenario_name)
