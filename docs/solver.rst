How the solver works and how to change equations
================================================

Changing soil, climate or rock inputs usually only requires editing
:doc:`the example notebook <_generated/Example>`. To change a chemical
relationship, start in ``smew/equations.py``.

How each timestep is solved
----------------------------

SMEW updates the stored amounts of water and chemical species from their
fluxes, then solves the coupled chemical-equilibrium equations. It starts
from an estimate based on the previous state and adjusts concentrations and
exchange fractions until the equations are satisfied closely enough.

Mathematically, it seeks :math:`F(x; q)=0`, where :math:`x` contains the
unknowns and :math:`q` the fixed parameters and inventories. Each component
of :math:`F` is a *residual*: the mismatch in one balance or equilibrium
relationship. In the code, ``p`` holds the trial unknowns :math:`x`, and
``out`` holds their residuals.

Python mode uses SciPy's ``fsolve``; compiled mode uses cminpack's ``hybrd``.
Both use a modified Powell hybrid method, estimate derivatives numerically,
and share the equations in ``equations.py``.

Convergence and accuracy
-------------------------

Solver convergence means the numerical stopping criterion was met; it does
not establish physical accuracy. Check the equation residuals, mass balances
and sensitivity to the timestep ``dt``.

.. note::

   The main chemistry solve uses ``xtol=1e-12`` as a stopping tolerance for
   relative changes between successive estimates. A residual above ``0.1`` in
   absolute value triggers
   a retry with a different hydrogen estimate and ``xtol=1e-14``; residuals
   still above ``0.1`` cause ``CHEMISTRY_NO_CONVERGENCE``. Residuals use the
   model's scaled units, so rescaling an equation changes this check.
   These tolerances do not set the time resolution.

Trial concentrations are not constrained to be positive, so an intermediate
estimate can make an equation undefined. See :doc:`errors` for troubleshooting.

Find the equation to change
---------------------------

Edit the residual functions in ``smew/equations.py``. For example, this line
calculate the 8th residual:

.. literalinclude:: ../smew/equations.py
   :language: python
   :start-at:     out[7] = Ca*nZrs1000
   :end-at:     out[7] = Ca*nZrs1000
   :dedent: 4

Check your change
------------------

Use an :doc:`editable installation <installation>` and start in
:doc:`Python mode <usage>`.

1. Edit the residual in ``smew/equations.py`` and, where possible, evaluate it
   at a known solution to check signs, units and variable order.
2. Restart the notebook kernel and rerun a familiar scenario with the same
   weather inputs. Compare mass balances, concentrations and exchange fractions.
3. Before rebuilding, declare any new intermediate variables in
   ``equations.pxd`` as needed for compilation. New parameters or unknowns
   require the additional changes in :doc:`reference/internals`.
4. Rebuild the C extensions, then restart in compiled mode and compare with
   the Python results, allowing for small rounding differences:

   .. code-block:: bash

      # Rebuild after changing the equations.
      make install-native

See :doc:`installation` for the command without Make. Run ``make test`` to
check the existing scenarios without replacing their reference results.
