Installation and builds
=======================

Install from PyPI
-----------------

As in the README, the basic installation is:

.. code-block:: console

   pip install smew

Use Python 3.11 or newer. If a ready-to-install package (called a *wheel*) is
available for your system and Python version, it already contains the fast
versions of the model equations and the cminpack solver used to solve them.
You do not need a C compiler to install it. Numba is installed automatically;
it translates numerical functions into fast machine code the first time you
use them. See :doc:`usage` to choose how the model runs and which solver it uses.

For an isolated environment, use Python's built-in ``venv``:

.. code-block:: bash

   # Create a separate Python environment for SMEW.
   python -m venv .venv

   # Activate the environment: choose the command for your terminal.
   source .venv/bin/activate          # Linux or macOS
   .\.venv\Scripts\Activate.ps1      # Windows PowerShell
   .venv\Scripts\activate.bat        # Windows Command Prompt

   # Install SMEW in the activated environment.
   pip install smew

If no ready-to-install package is available for your system and Python version,
pip tries to build SMEW on your computer. This requires a C compiler. Use the
command below to install only a ready-to-install package; if none is available,
pip stops and reports an error:

.. code-block:: console

   python -m pip install --only-binary=smew smew

Editable development installation
---------------------------------

Clone the repository and run these commands from its root:

.. code-block:: bash

   # Download the SMEW source code.
   git clone https://github.com/MatteoBertagni/SMEW.git

   # Open the downloaded project folder.
   cd SMEW

   # Follow the commands from above to create and activate your
   # virtual environment

   # Install SMEW so that edits to its Python files take effect,
   # along with tools for testing and building the documentation.
   pip install -e '.[dev,docs]'

The ``dev`` extra installs the testing, packaging and marimo tools; ``docs``
installs Sphinx and the Jupyter notebook execution dependencies. Omit either
extra when it is not needed. Building the C extensions requires a C compiler
(GCC/Clang on Linux or macOS, or MSVC Build Tools on Windows).
Pip installs Cython and the other build dependencies in
an isolated build environment. cminpack is bundled in the repository; no
separate MINPACK library is required.

.. note::

   On Linux, you may also need to install Python's development files to build
   the C extensions. On Ubuntu or Debian, these are provided by ``python3-dev``
   (or a version-specific package such as ``python3.11-dev``). Use the package
   that matches your Python version.

Ordinary Python edits become available in the editable installation. Restart
the interpreter or notebook kernel after edits, since imported modules and
Numba specializations remain in memory. Changes to the Cython-compiled
equations or native solver also need a rebuild:

.. code-block:: console

   make install-native

The equivalent command, including on platforms without Make, is:

.. code-block:: console

   pip install --no-cache-dir -e .

This uses the default ``SMEW_BUILD_NATIVE=1``. Restart Python after rebuilding
to load the new extensions. Edit ``smew/equations.py`` and its declarations in
``smew/equations.pxd``; generated files under ``build/`` are overwritten.

Python mode for debugging and experiments
--------------------------------------------

Use Python mode to follow a calculation step by step, inspect intermediate
values, or try changes to the model equations. It is slower than the compiled
mode, but makes it easier to understand and investigate the model's behaviour.

If SMEW is already installed, you can select Python mode without reinstalling
anything. Put the following in the **first cell of your notebook**, or at the
top of your script, before importing Numba or SMEW:

.. code-block:: python

   import os

   # Run the model as Python code, using SciPy to solve the equations.
   os.environ["NUMBA_DISABLE_JIT"] = "1"
   import smew

Restart the notebook kernel and run the cells from the top after changing
this setting. Also restart after editing SMEW's Python files so that your
changes are loaded. For a script, start a new run. See :doc:`usage` for more
about choosing the solver mode.

To edit the model itself, download the source code and activate your environment
as described above. You can install this local copy without a C compiler.
From the SMEW folder, run:

.. code-block:: bash

   # Install your local copy so that you can edit its Python files.
   make install-python

If Make is unavailable, use the command for your terminal instead.
On Linux or macOS:

.. code-block:: bash

   # Skip building C extensions and install your local copy for editing.
   SMEW_BUILD_NATIVE=0 pip install --no-cache-dir -e .

On Windows PowerShell:

.. code-block:: powershell

   # Skip building C extensions during this installation.
   $env:SMEW_BUILD_NATIVE = "0"
   pip install --no-cache-dir -e .

   # Clear the installation setting afterwards.
   Remove-Item Env:SMEW_BUILD_NATIVE

These installation commands skip the C extensions; they do not select Python
mode. Keep ``NUMBA_DISABLE_JIT`` set to ``"1"`` in your notebook or script as
shown above. Numba remains installed, but is disabled by this setting.

.. note::

   If you previously built C extensions in this folder, their files may remain.
   They do not need to be removed to use Python mode. To check an installation
   with no C extensions at all, run ``make clean-build`` before
   ``make install-python``. This removes generated extensions and package build
   files, including the ``build/`` and ``dist/`` folders.

Makefile commands
-----------------

.. note::

   If ``make`` is unavailable, install it with ``sudo apt install make`` on
   Ubuntu/Debian, or install Apple's `Command Line Tools
   <https://developer.apple.com/documentation/xcode/installing-the-command-line-tools/>`_
   with ``xcode-select --install`` on macOS.

   Make is optional. The direct PowerShell commands shown on this page offer
   an alternative on Windows, with no Make installation needed.

Activate your Python environment, then run Make from the repository root.
Make uses the ``python`` command from your active environment. To choose
another interpreter, pass its command or path, for example
``make build PYTHON=python3.11``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Purpose
   * - ``make build`` / ``make build-native``
     - Build a native source distribution and wheel under ``dist/``.
   * - ``make build-python``
     - Build a ready-to-install Python package without C extensions.
   * - ``make install-native``
     - Rebuild extensions in the editable installation.
   * - ``make install-python``
     - Reinstall the editable project without building extensions.
   * - ``make docs`` / ``make docs-check``
     - Build HTML and execute a copy of ``Example.ipynb``; ``docs-check``
       rereads all sources. See :doc:`contributing`.
   * - ``make test``
     - Run Python and compiled suites in separate processes, including
       full-year regressions.
   * - ``make test-serial``
     - Run the same suites without parallel Python workers.
   * - ``make example``
     - Open the existing interactive marimo example.
   * - ``make check-notebooks``
     - Check the marimo example's structure; does not execute ``Example.ipynb``.
   * - ``make update-baseline CASE=default``
     - Intentionally replace one regression baseline after reviewing a change.
   * - ``make clean-build``
     - Remove generated native extensions, ``build/``, ``dist/`` and egg metadata.
   * - ``make clean-docs``
     - Remove ``docs/_build/``, ``docs/_generated/`` and ``docs/.jupyter_cache/``.
   * - ``make clean``
     - Run both ``clean-build`` and ``clean-docs``.

For documentation and packaging without Make, activate your Python environment
and run the following commands from the SMEW folder. These work in Windows
PowerShell as well as Linux and macOS terminals; choose the command for your task:

.. code-block:: powershell

   # Build the documentation website (equivalent to make docs).
   python -m sphinx -b html -W --keep-going docs docs/_build/html

   # Rebuild all pages and rerun the example (equivalent to make docs-check).
   python -m sphinx -b html -W --keep-going -E -a -D nb_execution_mode=force docs docs/_build/html

   # Build the standard installation packages (equivalent to make build).
   python -m build
