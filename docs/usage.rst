Using the model and choosing a solver
=====================================

Start from :doc:`the example notebook <_generated/Example>` and change its
soil, climate or rock inputs for your scenario.

Prepare a simulation
---------------------

Use the same time grid for all inputs, with ``dt`` in days. The notebook
follows this order:

1. Set temperature and rainfall, then calculate potential evapotranspiration
   with :mod:`smew.hydroclimatic`.
2. Calculate vegetation with :func:`smew.vegetation.veg` and soil moisture,
   infiltration and water losses with :func:`smew.moisture.moisture_balance`.
3. Calculate soil organic carbon and respiration with
   :func:`smew.organic_carbon.respiration`.
4. Set the initial chemistry with :mod:`smew.ic`: use ``f_CEC_to_conc`` when
   exchange fractions and pH are known, or ``conc_to_f_CEC`` when dissolved
   concentrations and pH are known.
5. Run :func:`smew.biogeochem.biogeochem_balance` with the rock properties.
   The result contains time series such as ``result["pH"]``, ``result["Alk"]``
   and ``result["M_rock"]``.

Choose the mode
----------------

**Python mode** is useful for inspecting calculations and testing changes to
an equation. **Compiled mode** runs faster for longer simulations or many
scenarios. Both use the same equation source.

Put this before importing SMEW or Numba, in the first notebook cell or at the
top of your script:

.. code-block:: python

   import os

   # "1": Python mode for experiments; "0": compiled mode for speed.
   os.environ["NUMBA_DISABLE_JIT"] = "1"
   import smew

Restart the notebook kernel and run all cells after changing the setting;
for a script, start a new run. Without this setting, SMEW uses compiled mode.
Python mode also works with an installation without C extensions; compiled
mode requires them. See :doc:`installation` for installation options.

The first compiled calculation in each Python session includes preparation
time, so run it once before measuring execution speed.

Comparisons
-----------

When comparing modes or equation changes, use the **same temperature and
rainfall time series** and allow for small rounding differences. The modes
generate random rainfall separately, so the same random seed alone does not
ensure identical rain. Check mass balances and sensitivity to ``dt`` as well
as whether the solver converges.

For warnings or failed calculations, see :doc:`errors`.
