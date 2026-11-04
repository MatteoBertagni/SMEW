"""Portable development commands for the setuptools-based SMEW project."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import click


ROOT = Path(__file__).resolve().parents[1]
FORWARD_ARGS = {"ignore_unknown_options": True}


def _run(module, *args, env=None):
    """Use the active interpreter and keep environment changes in the child."""
    command = [sys.executable, "-m", module, *args]
    click.echo(f"$ {subprocess.list2cmdline(command)}")
    result = subprocess.run(command, cwd=ROOT, env={**os.environ, **(env or {})})
    if result.returncode:
        raise click.exceptions.Exit(result.returncode)


def _build(native, args):
    _run("build", *args, env={"SMEW_BUILD_NATIVE": "1" if native else "0"})


def _install(native):
    _run("pip", "install", "--no-cache-dir", "-e", ".",
         env={"SMEW_BUILD_NATIVE": "1" if native else "0"})


def _test(args, *, parallel):
    workers = (
        ["-n", "auto", "--maxprocesses=3", "--dist=loadgroup"]
        if parallel else ["-n", "0"]
    )
    _run("pytest", *workers, "tests", "--ignore=tests/test_backend.py", *args,
         env={"NUMBA_DISABLE_JIT": "1"})
    _run("pytest", "-n", "0", "tests", *args,
         env={"NUMBA_DISABLE_JIT": "0", "NUMBA_NRT_STATS": "1"})


def _docs(args, *, force):
    rebuild = ["-E", "-a", "-D", "nb_execution_mode=force"] if force else []
    _run("sphinx", "-b", "html", "-W", "--keep-going", *rebuild, *args,
         "docs", "docs/_build/html")


def _remove(paths):
    for path in paths:
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)


def _clean_build():
    _remove(ROOT / name for name in ("build", "dist", "smew.egg-info", "src/smew.egg-info"))
    native = ROOT / "src" / "smew" / "_native"
    for pattern in (
        "_equations*.so", "_equations*.pyd", "_minpack*.so", "_minpack*.pyd", "_minpack.c",
    ):
        _remove(native.glob(pattern))


def _clean_docs():
    _remove(ROOT / "docs" / name for name in ("_build", "_generated", ".jupyter_cache"))


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def build(args):
    """Build a native source distribution and wheel under dist/."""
    _build(True, args)


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def build_native(args):
    """Build packages with C extensions (same as build)."""
    _build(True, args)


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def build_python(args):
    """Build packages without C extensions."""
    _build(False, args)


@click.command()
def install():
    """Rebuild and install C extensions in the editable project."""
    _install(True)


@click.command()
def install_native():
    """Install the editable project with C extensions (same as install)."""
    _install(True)


@click.command()
def install_python():
    """Install the editable project without building C extensions."""
    _install(False)


@click.command()
def clean():
    """Remove generated package and documentation files."""
    _clean_build()
    _clean_docs()


@click.command()
def clean_build():
    """Remove generated packages, native extensions and build files."""
    _clean_build()


@click.command()
def clean_docs():
    """Remove generated documentation and cached notebook results."""
    _clean_docs()


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def test(args):
    """Run both model backends, with parallel Python workers."""
    _test(args, parallel=True)


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def test_serial(args):
    """Run both model backends without parallel workers."""
    _test(args, parallel=False)


@click.command()
def example():
    """Open the interactive marimo example."""
    _run("marimo", "edit", "tests/example_marimo_notebook.py")


@click.command()
def check_notebooks():
    """Check the marimo example's structure without running it."""
    _run("marimo", "check", "--strict", "tests/example_marimo_notebook.py")


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def docs(args):
    """Build the documentation website and run the example notebook."""
    _docs(args, force=False)


@click.command(context_settings=FORWARD_ARGS)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
def docs_check(args):
    """Rebuild all documentation and rerun the example notebook."""
    _docs(args, force=True)


@click.command()
@click.option("--case", required=True, help="Scenario name from tests/checks.toml.")
def update_baseline(case):
    """Replace one regression baseline after reviewing a model change."""
    _run("tests.update_baseline", "--case", case)
