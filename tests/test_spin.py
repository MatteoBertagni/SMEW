"""Check command portability without building or running model scenarios."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tomllib
from unittest.mock import Mock

from click.testing import CliRunner
import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("smew_spin", ROOT / ".spin" / "cmds.py")
commands = importlib.util.module_from_spec(spec)
spec.loader.exec_module(commands)


@pytest.fixture
def run(monkeypatch):
    mock = Mock(return_value=subprocess.CompletedProcess([], 0))
    monkeypatch.setattr(commands.subprocess, "run", mock)
    return mock


def test_registered_commands_load():
    with (ROOT / "pyproject.toml").open("rb") as file:
        config = tomllib.load(file)
    result = subprocess.run(
        [sys.executable, "-m", "spin", "--help"], cwd=ROOT,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "Could not" not in result.stdout
    for section in config["tool"]["spin"]["commands"].values():
        for entry in section:
            name = getattr(commands, entry.split(":")[1]).name
            assert name in result.stdout


@pytest.mark.parametrize("parallel", [True, False])
def test_backends_use_separate_environments_and_forward_arguments(run, monkeypatch, parallel):
    monkeypatch.setenv("NUMBA_DISABLE_JIT", "parent")
    monkeypatch.setenv("NUMBA_NRT_STATS", "parent")
    command = commands.test if parallel else commands.test_serial
    result = CliRunner().invoke(command, ["--", "-k", "physical or regression"])
    assert result.exit_code == 0, result.output
    python, compiled = run.call_args_list
    for call in (python, compiled):
        assert call.args[0][:3] == [sys.executable, "-m", "pytest"]
        assert call.args[0][-2:] == ["-k", "physical or regression"]
        assert call.kwargs["cwd"] == ROOT
        assert call.kwargs["env"]["PATH"] == os.environ["PATH"]
    assert python.kwargs["env"]["NUMBA_DISABLE_JIT"] == "1"
    assert compiled.kwargs["env"]["NUMBA_DISABLE_JIT"] == "0"
    assert compiled.kwargs["env"]["NUMBA_NRT_STATS"] == "1"
    assert "--ignore=tests/test_backend.py" in python.args[0]
    assert "--ignore=tests/test_backend.py" not in compiled.args[0]
    assert ("auto" in python.args[0]) == parallel
    assert compiled.args[0][3:5] == ["-n", "0"]
    assert os.environ["NUMBA_DISABLE_JIT"] == "parent"
    assert os.environ["NUMBA_NRT_STATS"] == "parent"


@pytest.mark.parametrize("failure_at", [0, 1])
def test_test_failure_is_propagated_and_stops_later_work(run, failure_at):
    run.side_effect = [subprocess.CompletedProcess([], 0)] * failure_at + [
        subprocess.CompletedProcess([], 7),
    ]
    result = CliRunner().invoke(commands.test)
    assert result.exit_code == 7
    assert run.call_count == failure_at + 1


@pytest.mark.parametrize("name,native", [
    ("build", "1"), ("build_native", "1"), ("build_python", "0"),
    ("install", "1"), ("install_native", "1"), ("install_python", "0"),
])
def test_build_mode_overrides_parent_without_changing_it(run, monkeypatch, name, native):
    monkeypatch.setenv("SMEW_BUILD_NATIVE", "parent")
    args = ["--", "--outdir", "folder with spaces"] if name.startswith("build") else []
    result = CliRunner().invoke(getattr(commands, name), args)
    assert result.exit_code == 0, result.output
    call = run.call_args
    assert call.kwargs["env"]["SMEW_BUILD_NATIVE"] == native
    assert call.args[0][:2] == [sys.executable, "-m"]
    if args:
        assert call.args[0][2:] == ["build", "--outdir", "folder with spaces"]
    else:
        assert call.args[0][2:] == ["pip", "install", "--no-cache-dir", "-e", "."]
    assert os.environ["SMEW_BUILD_NATIVE"] == "parent"


@pytest.mark.parametrize("force", [False, True])
def test_sphinx_arguments_and_forced_execution(run, force):
    command = commands.docs_check if force else commands.docs
    result = CliRunner().invoke(command, ["--", "-D", "language=fr"])
    assert result.exit_code == 0, result.output
    argv = run.call_args.args[0]
    assert argv[:7] == [sys.executable, "-m", "sphinx", "-b", "html", "-W", "--keep-going"]
    assert argv[-4:] == ["-D", "language=fr", "docs", "docs/_build/html"]
    assert ("nb_execution_mode=force" in argv) == force
    assert ("-E" in argv) == force
    assert ("-a" in argv) == force


def test_baseline_requires_explicit_case(run):
    result = CliRunner().invoke(commands.update_baseline)
    assert result.exit_code == 2
    run.assert_not_called()
    result = CliRunner().invoke(commands.update_baseline, ["--case", "default"])
    assert result.exit_code == 0
    assert run.call_args.args[0][2:] == ["tests.update_baseline", "--case", "default"]


@pytest.mark.parametrize("name,args", [
    ("example", ["edit", "tests/example_marimo_notebook.py"]),
    ("check_notebooks", ["check", "--strict", "tests/example_marimo_notebook.py"]),
])
def test_notebook_commands(run, name, args):
    result = CliRunner().invoke(getattr(commands, name))
    assert result.exit_code == 0
    assert run.call_args.args[0] == [sys.executable, "-m", "marimo", *args]


@pytest.mark.parametrize("name", ["clean", "clean_build", "clean_docs"])
def test_cleanup_is_scoped_and_repeatable(tmp_path, monkeypatch, name):
    monkeypatch.setattr(commands, "ROOT", tmp_path)
    build_files = [
        "build/object.o", "dist/package.whl", "smew.egg-info/PKG-INFO",
        "src/smew.egg-info/PKG-INFO", "src/smew/_native/_equations.test.so",
        "src/smew/_native/_equations.test.pyd", "src/smew/_native/_minpack.test.so",
        "src/smew/_native/_minpack.test.pyd", "src/smew/_native/_minpack.c",
    ]
    docs_files = ["docs/_build/index.html", "docs/_generated/Example.ipynb", "docs/.jupyter_cache/cache"]
    source_files = [
        "src/smew/equations.py", "src/smew/_native/_minpack.pyx",
        "src/smew/_native/__init__.py", "src/smew/_native/unrelated.pyd",
        "docs/installation.rst", "notebooks/Example.ipynb", "tests/baselines/default.npz",
    ]
    for filename in build_files + docs_files + source_files:
        path = tmp_path / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("keep or remove", encoding="utf-8")
    for _ in range(2):
        result = CliRunner().invoke(getattr(commands, name))
        assert result.exit_code == 0, result.output
    for filename in build_files:
        assert (tmp_path / filename).exists() == (name == "clean_docs")
    for filename in docs_files:
        assert (tmp_path / filename).exists() == (name == "clean_build")
    for filename in source_files:
        assert (tmp_path / filename).read_text(encoding="utf-8") == "keep or remove"
