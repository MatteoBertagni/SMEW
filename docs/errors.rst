When a simulation fails
=======================

SMEW reports a numerical failure with an error name, a message and, where
available, the timestep and conditions at which it occurred.

Read the message
-----------------

For example, ``CHEMISTRY_NO_CONVERGENCE`` means that the chemistry equations
were not satisfied closely enough, even after a second attempt.
``CHEMISTRY_RESIDUAL`` can indicate an undefined value during the calculation,
such as division by zero.

Look for ``time_days`` to locate the failure. ``previous_pH`` is the last
completed state's pH; ``trial_H`` is a value tried by the solver, not an accepted
concentration. ``nan`` means that a numerical value is undefined. The step
number starts at zero.

Check the simulation
---------------------

1. Check input units, array lengths and the order of chemical species against
   :doc:`usage`. Exchange fractions should sum to one.
2. Check the initial chemistry and inputs near the failed timestep for
   inconsistent or physically implausible values.
3. Try a smaller ``dt`` with the same physical forcing. Recalculate all inputs
   on the new time grid and check whether the results become stable as ``dt``
   decreases.
4. If you changed an equation, use Python mode to inspect its values and
   compare with a known case; see :doc:`solver`.

A solver warning does not always stop the simulation, but inspect the results
and equation balances before relying on them. Tightening the solver tolerance
alone does not demonstrate physical accuracy.

Printing errors
---------------

The example below displays the error name and message, then stops the current run:

.. code-block:: python

   import smew

   try:
       result = smew.biogeochem_balance(**inputs)  # Your model inputs.
   except smew.SMEWError as error:
       print(error.code.name, error.message)
       raise  # Stop this run; no completed result is available.

See :doc:`reference/errors` for the error-code reference and instructions for
calling SMEW from another compiled model.
