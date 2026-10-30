Solver implementation reference
================================

Use this page when adding a parameter, an unknown, or a new set of equations.
For changes to an existing formula and checks on the results, see
:doc:`../solver`.

Files to edit
--------------

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - Source file
     - What it contains
   * - `smew/equations.py <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/equations.py>`_
     - Equations, number of unknowns, and functions used by SciPy.
   * - `smew/equations.pxd <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/equations.pxd>`_
     - Variable and argument types needed to compile the equations.
   * - `smew/biogeochem.py <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/biogeochem.py>`_ / `smew/ic.py <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/ic.py>`_
     - Model inputs, initial estimates and use of the calculated solution.
   * - `smew/_utils.py <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/_utils.py>`_
     - Selection of the equation system in Python mode.
   * - `smew/_native/_minpack.pyx <https://github.com/MatteoBertagni/SMEW/blob/HEAD/smew/_native/_minpack.pyx>`_
     - The same selection in compiled mode, including parameter and unknown counts.

Changes to keep consistent
---------------------------

**Adding a parameter:** pass it in the same order through the model, equations
and both solvers. Update its type in ``equations.pxd`` and the parameter count
in ``_system_layout`` in ``_minpack.pyx``.

**Adding an unknown:** add its equation and initial estimate, update the number
of unknowns and resize the solver's working arrays. Keep the variable order
consistent when reading inputs and extracting results, with one equation per
unknown. Update the compiled declarations and ``_system_layout`` as needed.

**Adding an equation system:** define its equations, register it in ``_utils.py``
and ``_minpack.pyx``, and add its declarations in ``equations.pxd``. Functions
that create Python arrays may also need an entry in ``setup.py``'s
``python_only`` list. Test the new system in both Python and compiled modes.

Edit these source files, not the generated copies under ``build/``. Rebuild
and compare results as described in :doc:`../solver`.
