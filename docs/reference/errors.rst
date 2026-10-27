Error codes and exceptions
==========================

In a notebook or Python script, SMEW reports numerical failures automatically
as ``SMEWError``. Its ``code.name`` identifies the problem and its ``message``
gives the conditions at which it occurred. See :doc:`../errors` for troubleshooting.

For repeated simulations, save the error name and message with the scenario
inputs so that you can identify and investigate failed runs later.

Calling from your own compiled model
-------------------------------------

This section is only needed when calling SMEW inside another Numba function.
The functions ``biogeochem_balance_numba``, ``respiration_numba``,
``conc_to_f_CEC_numba``, ``total_to_f_CEC_and_conc_numba`` and ``Kelland_numba``
return ``(result, code, message)`` instead of raising an error directly.

Code ``0`` means success; a nonzero code means failure and ``result`` is ``None``.
Return the code and message to ordinary Python, then call
``smew.raise_for_error(code, message)`` there. Initialization functions can also
return a warning message with code ``0``; this should be reported separately.

Python reference
-----------------

.. automodule:: smew.errors
   :members:
   :undoc-members:
