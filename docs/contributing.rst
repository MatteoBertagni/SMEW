Maintaining the documentation
=============================

Edit the README for the home page, ``docs/*.rst`` for the guides, or the
explanatory text inside Python functions (their *docstrings*) for the API
reference. Edit ``notebooks/Example.ipynb`` to change the example.

Preview your changes
---------------------

From the SMEW folder, with your Python environment activated:

.. code-block:: bash

   # Install the documentation tools once.
   pip install -e ".[docs]"

   # Build the pages and run a copy of the example notebook.
   spin docs

Open ``docs/_build/html/index.html`` in your browser.
For text-only changes, ``spin docs -- -D nb_execution_mode=off`` builds a
faster preview by skipping notebook execution.
Notebook plots and calculated results are omitted.
These commands work on Windows, macOS and Linux.
See :doc:`installation` for the full command list.

Check before sharing
---------------------

.. code-block:: bash

   # Rebuild every page and rerun the example to catch errors.
   spin docs-check

Warnings, broken internal references and notebook execution errors must be
resolved before the check passes. The notebook uses the mode selected in its
first cell; its Python-mode setting is currently commented out.

.. note::

   The build runs a temporary copy of the notebook and leaves the original
   unchanged. ``spin clean-docs`` removes generated pages and cached results
   if you need a fresh build.
