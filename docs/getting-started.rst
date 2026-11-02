Getting started
===============

Run a one-year enhanced weathering simulation with the example notebook.

* :doc:`Installation and builds <installation>` — environment setup,
  development installation and Make commands.
* :doc:`Example: one year of enhanced weathering <_generated/Example>` —
  the full notebook with explanations, calculations and plots.

1. Install SMEW
-----------------

Use Python 3.11 or newer. In your terminal, with your Python environment activated:

.. code-block:: bash

   # Install the model.
   pip install smew

   # Install JupyterLab if you do not already have a notebook application.
   pip install jupyterlab

2. Open the example
---------------------

Download `Example.ipynb from GitHub
<https://github.com/MatteoBertagni/SMEW/blob/HEAD/Example.ipynb>`_ using
**Download raw file**, or use the copy in the repository's main folder.
From the folder containing the notebook, run:

.. code-block:: bash

   # Start JupyterLab, then double-click Example.ipynb in its file browser.
   jupyter lab

3. Run the simulation
-----------------------

For Python mode, remove the leading ``#`` from this line in the first code
cell, keeping it before ``import smew``:

.. code-block:: python

   os.environ["NUMBA_DISABLE_JIT"] = "1"

In the **Kernel** menu, choose **Restart Kernel and Run All Cells**.
The notebook calculates the soil response and plots moisture, pH, alkalinity
and rock dissolution over the year.

.. toctree::
   :hidden:
   :maxdepth: 1

   installation
   _generated/Example
