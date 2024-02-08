Public GEMS API
========================

The module :mod:`gems` is organised in several submodules.

* :mod:`gems.libio` contains input/output routines. You want to start with this one
    in order to load the data and pass it to the fitting routines.
* :mod:`gems.fit` contains the optimisation routines for fitting the data.
* :mod:`gems.libuk` contains all the tools you need to handle the microspeciation
    information.
* :mod:`gems.plotter` contains routines for plotting the results obtained
* :mod:`gems.report` contains reporting routines

variable definition
************************

.. _def_free_energy:

Free Energy
------------------------
:py:class:`dict` where the keys are the microstate codified as a tuple of ones and zeros
and the value is a float representing the value of the energy that correspond to 
that microstate.

.. code-block:: python

    free_energy = {(0, 0, 0): 0.0, (0, 0, 1): -12.41, (0, 1, 0): -13.01,
        (0, 1, 1): -22.14, (1, 0, 0): -10.43, (1, 0, 1): -29.13,
        (1, 1, 0): -27.34, (1, 1, 1): -41.42}

.. _microstate_definition:

Microstate
----------
The microstate can be codified either with numbers or letters. A tuple of zeros and
ones of the lenght of the number of protonation centres can define each one of the
microstates. When the centre s:sub:`j` is protonated, a 1 occupies that place in the
tuple, zero otherwise. A given microstate is generally noted as {*s*:sub:`j`} which
means the collection of 0s and 1s that represent it, e.g. {0, 1, 0, 1} meaning that
*s*:sub:`1` = 0, *s*:sub:`2` = 1, *s*:sub:`3` = 0, *s*:sub:`4` = 1. 

.. _def_microstate_probability:

Microstate Probability
----------------------

Each microstate has a probability that can be calculated from the free energy.

.. math::
    p(\{s_j\}) = \frac{a_{\mathrm{H}}^n e^{-\beta F(\{s_j\})}}{\sum_{\{s_j\}} a_{\mathrm{H}}^n e^{-\beta F(\{s_j\})}}

Internally, this variable is stored as a :class:`dict` where the key is a
microstate binary tuple and the value is an iterable (like a 1D array) containing
the values of the microstate probability for each titration point.

.. _def_macrostate_probability:

Macrostate Probability
----------------------

Each macrostate has a probability that can be calculated according to

.. math::
    P_n(a_H) = \frac{\bar{K_n} a_{\rm H}^n}{\sum_{n=0}^N \bar{K_n} a_{\mathrm{H}}^n}


.. _def_conditional_probability:

Conditional Probability
-----------------------

.. math::
    \pi_n(\{s_j\}) = \frac{e^{-\beta F(\{s_j\})}}{\sum_{\{s_j\}} e^{-\beta F(\{s_j\})}}

.. math::
    \pi_n(\{s_j\}) = \pi_n({s_j}) P_n(a_H)


.. _def_molecule_symmetry:

Molecule Symmetry
-----------------
The symmetry of a molecule can be expressed as a string of letters with a multiplier for
the sites that are equivalent. For example :samp:`A2B` if the molecule has three
centres, two of which are equivalent; or :samp:`ABC` if the molecule has three 
different centres.

.. _def_a_matrix:

A Matrix
--------

The **A** matrix represents the variation in population (rows) when the


.. _def_b_matrix:

B Matrix
--------

The **B** matrix represents the variation of the chemical shift (rows) when the
macroscopic protonation degree changes in one unit.

.. math::
    \delta_l - \delta^{(0)}_l = \sum_{n=1}^N B_{ln} P_n(a_H)

.. _def_delta_matrix:
   
Delta Matrix
------------

The matrix Δ represents the variation of the chemical shift (rows) when the
protonation degree of a particular site (columns) goes from protonated to
unprotonated.

.. math::
    \delta_l - \delta^{(0)}_l = \sum_{m=1} \Delta_{lm} \theta_m


libuk.py
********

.. automodule:: gems.libuk
    :members:
    :undoc-members:
    :show-inheritance:

report.py
*********

.. automodule:: gems.report
    :members:
    :undoc-members:
    :show-inheritance:

fit.py
******

.. automodule:: gems.fit
    :members:
    :undoc-members:
    :show-inheritance:


plotter.py
**********

.. automodule:: gems.plotter
    :members:
    :undoc-members:
    :show-inheritance:

libio.py
********

.. automodule:: gems.libio
    :members:
    :undoc-members:
    :show-inheritance:
