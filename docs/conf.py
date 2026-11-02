"""Sphinx configuration; import the Python API without compiling it for autodoc."""

import os
from pathlib import Path
import sys

from docutils import nodes
from sphinx import addnodes
import nbformat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Plain functions expose signatures/docstrings to autodoc, including public
# aliases of @njit functions. Restore the environment before notebook kernels
# start; each example selects its own backend in its first code cell.
_jit_setting = os.environ.get("NUMBA_DISABLE_JIT")
os.environ["NUMBA_DISABLE_JIT"] = "1"
import smew
if _jit_setting is None:
    os.environ.pop("NUMBA_DISABLE_JIT", None)
else:
    os.environ["NUMBA_DISABLE_JIT"] = _jit_setting

project = "SMEW"
author = "SMEW contributors"
release = smew.__version__
version = release
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_nb",
]
root_doc = "index"
exclude_patterns = ["_build", ".jupyter_cache", "**/.ipynb_checkpoints", "diagrams"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": True}
napoleon_numpy_docstring = True
napoleon_google_docstring = True
myst_enable_extensions = ["html_image"]
myst_heading_anchors = 3

# Always execute: a cached notebook would not notice a change to the installed
# SMEW implementation. Errors fail the build.
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
nb_execution_in_temp = True

html_theme = "pydata_sphinx_theme"
html_title = f"SMEW {release} documentation"
html_logo = "../images/SMEW_logo.png"
html_theme_options = {
    "logo": {"alt_text": "SMEW documentation home"},
    "github_url": "https://github.com/MatteoBertagni/SMEW",
    "navbar_align": "left",
    "show_nav_level": 1,
    "navigation_depth": 3,
    "show_toc_level": 2,
}
html_static_path = ["_static"]
html_css_files = ["custom.css"]


def _readme_links(app, doctree):
    """Resolve repository-relative README links when included at the docs root."""
    if app.env.docname != "index":
        return
    for node in doctree.findall(addnodes.pending_xref):
        target = node.get("reftarget", "")
        if target == "Example.ipynb":
            node["reftarget"] = "_generated/Example"
            node["refdomain"] = "std"
            node["reftype"] = "doc"
        elif target.startswith("docs/"):
            node["reftarget"] = str(Path(target.removeprefix("docs/")).with_suffix(""))
            node["refdomain"] = "std"
            node["reftype"] = "doc"
        elif target in {"tests/README.md", "third_party/cminpack/CopyrightMINPACK.txt"}:
            link = nodes.reference(
                "", "", refuri=f"https://github.com/MatteoBertagni/SMEW/blob/HEAD/{target}"
            )
            link.extend(node.children)
            node.replace_self(link)
    # The README's centered HTML logo is copied locally by html_logo above.
    for node in doctree.findall(nodes.raw):
        if 'src="images/SMEW_logo.png"' in node.astext():
            node[0] = nodes.Text(node.astext().replace(
                'src="images/SMEW_logo.png"', 'src="_static/SMEW_logo.png"'
            ))


def _prepare_example(app):
    """Stage the root notebook without modifying its source or saved outputs."""
    notebook = nbformat.read(ROOT / "Example.ipynb", as_version=4)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    target = Path(app.srcdir) / "_generated" / "Example.ipynb"
    target.parent.mkdir(exist_ok=True)
    content = nbformat.writes(notebook)
    if not target.exists() or target.read_text(encoding="utf-8") != content:
        target.write_text(content, encoding="utf-8")


def setup(app):
    app.connect("builder-inited", _prepare_example)
    app.connect("doctree-read", _readme_links)
