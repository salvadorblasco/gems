Algorithm Details
=================

Microstate Nomenclature
***********************

A molecule with *N* protonation sites can have *H* protons
arranged in :math:`N \choose H` ways. The protonation
centre *j* can be either protonated *s*:sub:`j` =1 or
deprotonated *s*:sub:`j` =0. Therefore, a microstate can be 
defined by a collection of binary parameters {*s*:sub:`j`} such as
{0, 1, 0, 0}.

.. seealso:: :ref:`microstate_definition`.


Cluster Expansion
*****************

The algorithm is largely based on the Cluster Expansion Method by Borkovec and 
Koper ([Borkovec2000]_ and [Borkovec2002]_) that also incorporates the simplifications
arising from the molecule symmetry from Szakács and Noszál ([Szakacs1999]_).
In this approach, the free energy *F* can be decomposed
into the sum of parameters of first, second and third order.

First-order parameters can be assimilated to the *raw* protonation constant of
a particular site neglecting any external influence.

.. image:: media/params1.png
   :scale: 50 %
   :alt: First-order parameters

Second-order parameters quantify the influence in the change of the microconstant
of a particular site when the neighbour site is protonated but the rest of
neighbouring sites are not. We assume that the 
influence is reciprocal and equal: :math:`\varepsilon_{ij} = \varepsilon_{ji}`
and that the self-influcence is zero: :math:`\varepsilon_{ii} = 0`.

.. image:: media/params2.png
   :scale: 50 %
   :alt: Second-order parameters

Third-order parameters quantify the change in second-order parameters when the
neighbour of a neighbour is protonated. We assume again that the reciprocal
influence is the same: 
:math:`\lambda_{ijk} = \lambda_{jik} = \lambda_{ikj} = \lambda_{kij} = \lambda_{jki} = \lambda_{kij}`
and the self-influence is zero: :math:`\lambda_{iii} = 0`, :math:`\lambda_{iij} = 0`, etc.

.. image:: media/params3.png
   :scale: 50 %
   :alt: Third-order parameters

The expansion could go on with higher-order parameters, however the absolute
value of the parameters approach zero as the distance increases and for 
a distance larger than two it can safely be considered zero.

The total number of parameters would be

* *n* first-order parameters
* :math:`{n \choose 2} = \frac{n(n-1)}2` Second-order parameters
* :math:`{n \choose 3} = \frac{n(n-1)(n-2)}3` Third-order parameters

As the system grows larger and larger, the number of microconstants grows
exponentially as :math:`2n^{(n-1)}`, but with the cluster expansion
capped to third-order only experiments a polynomial growth.
 
.. image:: media/params_growth.png
   :scale: 50 %
   :alt: Growth of the number of parameters


Symmetry Simplification
***********************

Most molecules have some kind of symmetry. This means that some
protonation centres are equivalent and, consequently, the number of
parameters is reduced because they can be constrained to have the
same value.

.. seealso:: :ref:`info_symmetry`


+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
| *n* | symmetry | without cluster expansion | with cluster expansion | with symmetry simplification | example                                                                        |
+=====+==========+===========================+========================+==============================+================================================================================+
|  2  | A2       |           3               |                 3      |             2                | ethylenediamine                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  2  | AB       |           3               |                 3      |             3                | alanine                                                                        |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  3  | A3       |          12               |                 7      |             3                | trisaminoethylamine                                                            |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  3  | A2B      |          12               |                 7      |             5                | citrate                                                                        |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  3  | ABC      |          12               |                 7      |             7                | inositol phosphate                                                             |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  4  | A4       |          32               |                14      |             3                | EDTA (carboxylate)                                                             |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  4  | A3B      |          32               |                14      |             6                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  4  | A2B2     |          32               |                14      |             7                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  4  | A2BC     |          32               |                14      |            10                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  4  | ABCD     |          32               |                14      |            14                | trilysine (amino groups)                                                       |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  5  | A5       |          80               |                25      |             3                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  5  | A4B      |          80               |                25      |             6                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  5  | ABCDE    |          80               |                25      |            25                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  6  | A6       |         192               |                41      |             3                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  6  | A5B      |         192               |                41      |             6                |                                                                                |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  6  | A4B2     |         192               |                41      |             8                | triethylenetetraminehexa-acetate (carboxylate)                                 |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 
|  6  | ABCDEF   |         192               |                41      |            41                |  corticotropin                                                                 |
+-----+----------+---------------------------+------------------------+------------------------------+--------------------------------------------------------------------------------+ 


Theoretical Background
**********************

With this decomposition we can fit the free energy of the system with
the following equation:

.. math:: 
    \frac{\beta F(\{s_j\})}{\ln 10} = -\sum_j{p\hat{K}_js_j} 
       + \frac1{2!}\sum_{ij}\varepsilon_{ij}s_is_j
       + \frac1{3!}\sum_{ijk}\lambda_{ijk}s_is_js_k

where *β* is inverse of the the thermal energy *k*:sub:`B`\ *T*.

It is possible to choose whether to fit only first order, both first and second
order or all of them. It is advisable to start with first order only and then
progress towards fitting the rest.  

With the free energy calculated, the macroconstants can be
derived from the following equation.

.. math::
    \bar{K_n} = \sum_{\{s_j\}} e^{-\beta F(\{s_j\})} \delta_{n,\sum_js_j}

And then, the probability of each macrostate is calculated

.. math::
    P_n(a_H) = \frac{\bar{K_n} a_{\rm H}^n}{\sum_{n=0}^N \bar{K_n} a_{\mathrm{H}}^n}

Then, experimental data is fitted to the equation

.. math::
    \delta_i = \delta_i^{(0)} + \sum_{m=1}^N B_{ln} P_n(a_{\mathrm{H}})

and the residual is calculated. Parameters are readjusted 
and the sequence starts over until the sum of the residuals squared are 
minimized.

.. math::
    p(\{s_j\}) = \Xi^{-1} a_{\mathrm{H}}^n e^{-\beta F(\{s_j\})} 

where

.. math::
    \Xi = \sum_{\{s_j\}} a_{\mathrm{H}}^n e^{-\beta F(\{s_j\})} 
        = \sum_{n=0}^N \bar{K_n} a_{\mathrm{H}}^n


.. math::
    \delta_i = \delta_i^{(0)} + \sum_{m=1}^N \Delta_{lm}\theta_m

.. math::
    \theta_m = \frac1N \sum_{\{s_j\}} s_m p(\{s_j\})
             = \sum_{n=0}^N A_{mn} P_n(a_{\mathrm{H}})

where

.. math::
   A_{mn} = \sum_{\{s_j\}} s_m \pi_{n}(\{s_j\}) \delta_{n, \sum_js_j}

The *microstate probability* can be defined as 

.. math::
    p(\{s_j\}) = \pi_n(\{s_j\}) P_n(a_{\rm H})
               = \Xi^{-1} a_H^n e^{-\beta F(\{s_j\})}

and the *conditional probability* is 

.. math::
   \pi(\{s_i\}) = \bar{K}_n^{-1} e^{-\beta F(\{s_j\})}



Algorithm Implementation
************************

The calculations are carried by SciPy in the background through the
function :func:`scipy.optimize.minimize`. The Minimisation
algorithm is the default one for the parameters used: the
L-BFGS-B algorithm. In the future, other algorithms will be available.
The objective function (see :func:`fit.fobj`) is passed and minimised.

References
**********

.. [Borkovec2000] *A Cluster Expansion Method for the Complete Resolution of Microscopic Ionization Equilibria
    from NMR Titrations*, Michal Borkovec and Ger J. M. Koper, *Anal. Chem.* **2000**\ , 72, 3272-3279.

.. [Borkovec2002] *Resolution of Microscopic Protonation Mechanisms in Polyprotic Molecules*, 
    Michal Borkovec, Marcin Brynda, Ger J. M. Koper and Bernard Spiess, *Chimia* **2002**, 56 695–701.

.. [Szakacs1999] *Protonation microequilibrium treatment of polybasic compounds with any possible symmetry*,
    Zoltán Szakács and Béla Noszál, *J. Math. Chem.* **1999**, 26, 139–155.
