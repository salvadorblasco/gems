# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                   !
# !                                                                          !
# ! Copyright (c) 2024 by Salvador Blasco <salvador.blasco@gmail.com>        !
# ! Licensed under MIT license (see file LICENSE)                            !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""This module contains the basic tools for microconstants. These functions
are loaded by :program:`GEMS` or they can be used in a standalone
python script.

Microstates Nomenclature
------------------------
* :func:`name_microstate`
* :func:`name_terms`
* :func:`expand_name`
* :func:`compact_name`
* :func:`molecule_info`

Functions for Multiplicity and Combinatory
------------------------------------------
* :func:`order_terms`
* :func:`kdeflation`
* :func:`comb`
* :func:`remainder`
* :func:`generate_microstep`
* :func:`state_microsteps`
* :func:`total_microconsts`
* :func:`total_microstates`
* :func:`generate_microstate`
* :func:`generate_all_microstates`
* :func:`count_microstate`
* :func:`remove_equivalent_states`

Terms calculations
------------------
* :func:`fit_free_energy`
* :func:`microstate_probability`
* :func:`macrostate_probability`
* :func:`matrix_a`
* :func:`matrix_b`
* :func:`conditional_probability2`
* :func:`conditional_probability1`
* :func:`avg_prtdeg`
* :func:`microstate_multiplicity`
* :func:`micro_constants2`
* :func:`micro_constants`
* :func:`error_macro_constants`
* :func:`macro_constants`
* :func:`num_prot_centres`
"""

import itertools
from itertools import product as iprod
import collections
import functools
import math
from math import exp
import re
import typing

import numpy as np
from numpy.linalg import lstsq
import scipy.optimize

import gems.isomorph
import gems.cmatrix


# ------------------------------------
#   Molecule and Symmetry Operations
# ------------------------------------

def compute_mapping(molecule, isomorphisms):
    size = len(molecule)
    mapping = {}
    for level in range(1+size):
        microstates = generate_microstates(size, level)
        my_set = gems.isomorph.clasify_microstates(microstates, isomorphisms)
        mapping.update({p: q for q in my_set for p in q})
    return mapping


def compute_connectivity(molecule):
    return gems.cmatrix.connectivity_matrix(molecule)


def compute_isomorphisms(molecule, connectivity):
    return gems.isomorph.find_isomorphisms(connectivity, molecule)


# ----------------------------
#   Microstates Nomenclature
# ----------------------------


def compact_name(input_name):
    """Compact a name.

    Parameters:
        input_name (str): an expanded molecule name
    Returns:
        str: the compact name with multipliers
    Example:
        >>> compact_name('AAABCC')
        A3BC2

    See Also:
        :func:`expand_name`

    """
    c = collections.Counter(input_name)
    return "".join(l + str(n) if n > 1 else l for l, n in c.items())


def expand_name(input_name):
    """Expand names with multipliers.

    When a site has high multiplicity it is more convenient to express it with
    a number rather than with repetition (eg. :samp:`A3B2` instead of
    :samp:`AAABB`. This routine converts the former in the latter.

    Parameters:
        input_name (str): a compact name
    Return:
        str: an expanded name
    See Also:
        :func:`compact_name`

    >>> expand_name('A3B3C')
    AAABBBC
    >>> expand_name('ABC')
    ABC
    >>> expand_name('A2BAAC')
    AABAAC
    """
    return re.sub(r'(\w)(\d+)', lambda x: x.group(1)*int(x.group(2)), input_name)


def name_equivalent_microstates(molecule: str, mapping: dict) -> dict:
    tt = []
    for v in set(mapping.values()):
        n = "".join(sorted(l for l, i in zip(molecule, v[0]) if i))
        m = min(sum(a*2**q for q, a in enumerate(b)) for b in v)
        tt.append((n, m, v))

    a,b,c = zip(*sorted(tt))
    retv = {}

    co = collections.Counter(a)
    q = itertools.chain(*[('*'*q for q in range(co[k])) for k in sorted(co)])
    for i, l in enumerate(a):
        retv[c[i]] = l + next(q)
    return retv


def reaction(ums1: tuple, ums2: tuple, molecule: str, msnames: dict):
    # breakpoint()
    react = msnames[ums1]
    prodt = msnames[ums2]
    proton = list(prodt)
    for l in react:
        proton.remove(l)
    return react, "".join(proton), prodt


# ------------------------------------------------
#   Functions for Multiplicity and Combinatory
# ------------------------------------------------


def compute_multiplicity(mapping):
    return {k: len(k) for k in set(mapping.values())}


def filter_by_macro(mapping, n):
    yield from (i for i in set(mapping.values()) if sum(i[0]) == n)


def generate_all_microstates(size):
    yield from iprod((0, 1), repeat=size)


def generate_microstates(m, level):
    """
    >>> generate_microstates(3, 1)
    {(1, 0, 0), (0, 0, 1), (0, 1, 0)}
    >>> generate_microstates(3, 2)
    {(1, 0, 1), (1, 1, 0), (0, 1, 1)}
    """
    return set(itertools.permutations(level*(1,)+(m-level)*(0,), m)) 


def comb(n, k):
    """Combinatorial for n elements  and k.

    Parameters:
        n (int): the number of elements
        k (int): the k combination, where :math:`n >= k`
    Returns:
        int: the combinatorial of n and k

    """
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))


def kdeflation(n):
    """Statistical distribution of equilibrium constants.

    Parameters:
        n (int): the number of protonations
    Yields:
        float: the deflation of the constants

    """
    def ratiocomb(n, k):
        return comb(n, k) / comb(n, k-1)
    for i in range(1, n+1):
        yield ratiocomb(n, i) / ratiocomb(n, 1)


# ----------------------
#   Terms calculations
# ----------------------


def avg_prtdeg(msprob):
    r"""Calculate the average protonation degree.

    This routine uses eq (6) from Borkovec et al. [Borkovec2000]_

    .. math::
        \theta_m = \frac1N \sum_{\{s_j\}} s_m p(\{s_j\})

    Parameters:
        msprob (dict): the :ref:`microstate probability <def_microstate_probability>`.
    Returns:
        numpy.ndarray: average protonation degree of each centre.

    """
    n = len(next(iter(msprob.keys())))
    return np.vstack(tuple(sum(ms[m] * p for ms, p in msprob.items())/n
                           for m in range(n)))


def conditional_probability1(micro_prob, macro_prob):
    r"""Calculate the :ref:`def_conditional_probability`.

    This routine uses eq (7) from Borkovec et al. [Borkovec2000]_

    .. math::
        \pi_n(\{s_j\}) = \pi_n({s_j}) P_n(a_H)

    Parameters:
        micro_prob (dict): The :ref:`def_microstate_probability`.
        macro_prob (:class:`numpy.ndarray`): The
            :ref:`def_macrostate_probability`.

    """
    return {s: mipr/macro_prob[:, sum(s)] for (s, mipr) in micro_prob.items()}


def conditional_probability2(free_enrg):
    r"""Calculate the :ref:`def_conditional_probability` of each :ref:`microstate_definition`.

    .. math::
        \pi_n(\{s_j\}) = \frac{e^{-\beta F(\{s_j\})}}
                              {\sum_{\{s_j\}}e^{-\beta F(\{s_j\})}
        \delta_{n,\sum_js_j}}

    """
    norm = macro_constants(free_enrg)
    return {s: exp(-energy)/norm[sum(s)] for (s, energy) in free_enrg.items()}


def macro_constants(free_enrg: dict):
    r"""Calculate the macro constants from the energy data.

    Applies equation below to calculate the macro constants.

    .. math::
        \bar{K_n} = \sum_{\{s_j\}} e^{-\beta F(\{s_j\})} \delta_{n,\sum_js_j}

    Parameters:
        free_enrg (dict): dict containing the free energy for each microstate.
        n_protcentr (int): The number of protonation centres

    """
    n_protcentr = num_prot_centres(free_enrg)
    retv = (sum(exp(-energy) for s, energy in free_enrg.items() if sum(s) == n)
            for n in range(n_protcentr+1))
    return tuple(retv)


def microsteps(size, mapping):
    msteps = {}

    for n in range(size):
        steps = set()
        ms1 = generate_microstates(size, n)
        ms2 = generate_microstates(size, 1+n)
        for a, b in itertools.product(ms1, ms2):
            if (r:=sum(abs(j-i) for i, j in zip(a,b))) == 1:
                if mapping[a] in msteps:
                    msteps[mapping[a]].add(mapping[b])
                else:
                    msteps[mapping[a]] = set([mapping[b]])
    return msteps


def error_macro_constants(free_enrg, error_freenrgy):
    """Calculate error in macroconstants.

    Parameters:
        free_enrg (dict): dict containing the :ref:`def_free_energy` for each microstate.
        error_free_enrgy (dict): dict containing the free energy for each
            microstate.

    .. seealso: :func:`macro_constants`

    """
    # eq (9)
    n_protcentr = num_prot_centres(free_enrg)
    return [sum(exp(-2*energy) * error_freenrgy[s]
                for s, energy in free_enrg.items()
                if sum(s) == n)
            for n in range(n_protcentr+1)]


def matrix_a(cond_prob, n_centres):
    r"""Compute matrix A.

    Compute :ref:`def_a_matrix` using the equation

    .. math::

        A_{mn} = \sum_{\{s_j\}} s_m \pi_n(\{s_j\}) \delta_{n,\sum_js_j}

    Parameters:
        cond_prob (dict): The :ref:`def_conditional_probability` for each
            microstate as output by :func:`conditional_probability1` or
            :func:`conditional_probability2`.
        n_centres (int): The number of protonation centres
        n_indepctres (int): The number of independent protonation centres.

    Returns:
        :class:`numpy.ndarray` : The A matrix as defined above.

    """
    mtra = np.zeros((n_centres, n_centres+1))
    for m, n in itertools.product(range(n_centres), range(n_centres+1)):
        mtra[m, n] = sum(s[m]*mipr for s, mipr in cond_prob.items() if sum(s) == n)
    return mtra


def matrix_b(macro_prob, shifts):
    r"""Least square solving for the equation below.

    Parameters:
        macro_prob (:class:`numpy.ndarray`):
        shifts (:class:`numpy.ndarray`): The experimental chemical shifts.

    .. math::
        \delta_i = \delta_i^{(0)} + \sum_{m=1}^N B_{ln} P_n(a_{\mathrm{H}})

    .. warning::
        As of now, NumPy does not support linear algebra for masked
        arrays. Therefore, if a masked array is passed, the linear algebra
        routines are not used. Instead, a least square minimisation is
        done. The result should be approximately the same and, for most
        cases the used will not perceive the difference.

    .. versionadded:: 0.5

    """
    improb = np.linalg.pinv(macro_prob)
    # bmatrix = improb @ shifts
    bmatrix = np.dot(improb, shifts)
    # if np.ma.is_masked(shifts):
    #     bmatrix = np.ma.dot(improb, shifts)
    #     # bmatrix = np.ma.dot(improb, shifts).data
    #     # bmatrix = mlstsq(macro_prob, shifts).data
    # else:
    #     bmatrix = np.dot(improb, shifts)
    #     # bmatrix = improb @ shifts             # not implemented for numpy 1.21
    return bmatrix


def micro_constant(microstep1, microstep2, free_energy):
    """Calculate the micro-equilibrium constant for two microstates.

    Given the free energy between two states, calculate
    .. math:: \\exp(-\\beta(F_2-F_1))
    """
    f1 = free_energy[microstep1]
    f2 = free_energy[microstep2]
    return -(f2-f1)/2.3025851


def micro_constants(msteps, free_enrg):
    return {ms1: {ms2: micro_constant(ms1, ms2, free_enrg) for ms2 in steps}
            for ms1, steps in msteps.items()}


# def micro_constants(molecule, free_enrg):
#     r"""Calculate the micro constants from the energy data.
# 
#     Parameters:
#         molecule (str): The molecule symmetry
#         free_enrg (dict): dict containing the :ref:`def_free_energy` for each microstate.
#     Returns:
#         dict: the values of the microconstants for each microstate and
#             microstep.
# 
#     """
#     retv = {}
#     for ms, energy in free_enrg.items():
#         retv2 = {}
#         for pos, bit in enumerate(ms):
#             if bit == 1:
#                 continue
#             key = molecule[pos]
#             if key in retv2:
#                 continue        # not covered
# 
#             newms = list(ms)
#             newms[pos] = 1
#             energy2 = free_enrg[tuple(newms)]
#             k = -(energy2-energy)/2.3025851
#             retv2[key] = k
#         retv[ms] = retv2
#     return retv


def macrostate_probability(macro_consts: np.ndarray, proton_activity: np.ndarray) -> np.ndarray:
    r"""Calculate probability of a macrostate.

    Apply the equations (8) and (10) from [Borkovec2000]_ in order to calculate
    the probability of the macrostates.

    .. math::
        P_n(a_H) = \frac{\bar{K_n} a_{\rm H}^n}{\sum_{n=0}^N
                   \bar{K_n} a_{\mathrm{H}}^n}

    Parameters:
        macro_consts (:class:`numpy.ndarray`): An array containing the values
            if the macroconstants.
        ah (:class:`numpy.ndarray`): The activity of protons.

    Returns:
        :class:`numpy.ndarray`: The probability of the macrostate.

    """
    n = np.arange(len(macro_consts))
    p = macro_consts[None, :]*proton_activity[:, None]**n[None, :]
    norm = np.sum(p, axis=1)
    return p / norm[:, None]


def merge_equivalent_microstates(mapping: dict, property_: dict, function: typing.Callable) -> dict:
    return {ums: functools.reduce(function, (property_[ms] for ms in ums))
            for ums in set(mapping.values())}


# DEPRECATED
# def microstate_population(conditional_probability, molecule):
#     """Calculate the population for each microstate.
# 
#     The population for each microstate is basically computed in the
#     :ref:`def_conditional_probability` variable. However, some microstates
#     are equivalent. This routine sums all equivalent microstates together and
#     returns a more sensible result.
# 
#     Parameters:
#         conditional_probability (dict): the :ref:`def_conditional_probability`.
#         molecule (str): A string representing the :ref:`def_molecule_symmetry`.
#     Returns:
#         dict: the microstate population.
# 
#     """
#     chain_ = itertools.chain.from_iterable
#     names = set(chain_(name_terms(molecule, i) for i in range(len(molecule)+1)))
# 
#     pop = {name: sum(conditional_probability[ms]
#                      for ms in conditional_probability
#                      if name_microstate(molecule, ms) == name)
#            for name in names}
#     return pop


def microstate_probability(free_energy: dict, ah: np.ndarray) -> dict:
    r"""Calculate microstate probability from the energy data.

    The probability of a given microstate is calculates from the free energy
    terms according to:

    .. math::
        p(\{s_j\}) = \Xi^{-1} a_H^n e^{-\beta F(\{s_j\})}

    Where :math:`\Xi = \sum_{\{s_j\}} a_H^n e^{-\beta F(\{s_j\})}`
    and :math:`n=\sum_js_j`

    Parameters:
        free_energy (dict): A dict containing the free energy for each
            microstate.
        ah (:class:`numpy.ndarray`): The activity of protons.

    Returns:
        :py:class:`dict`: The microstate probability for each microstate

    """
    p = {s: ah**sum(s) * math.exp(-energ) for s, energ in free_energy.items()}
    norm = sum(p.values())
    return {s: v/norm for s, v in p.items()}


def mlstsq(A, B):
    """Linear least squares for masked arrays.

    Compute the array *x* so that :math:`||Ax - B||` is mimized.

    Parameters:
        A (:class:`numpy.ndarray`):
        B (:class:`numpy.ma.MaskedArray`):

    Returns:
        :class:`numpy.ndarray` : The result

    .. seealso:: :func:`scipy.linalg.lstsq`
    .. warning:: This function will be removed when NumPy finally supports
        masked matrix inversion.

    """
    mres = lstsq(A, B, rcond=None)[0]
    for c in (n for n, t in enumerate(np.any(B.mask, axis=0)) if t):
        _mask = B.mask[:, c]
        b = B[:, c]
        retv = scipy.optimize.lsq_linear(A[~_mask, :], b[~_mask])
        mres[:, c] = retv['x']
    return mres


def num_prot_centres(data: dict) -> int:
    """Compute the number of protonation centres based on dict of data.

    Parameters:
        data (:class:`dict`): any type of data where the keys are a tuple
            of zeros and nones indicating the microstate.
    Returns:
        int: the number of protonation centres.

    """
    return len(tuple(data.keys())[0])


def num_prot_centres2(data: dict) -> int:
    """Compute the number of protonation centres based on dict of data.

    Parameters:
        data (:class:`dict`): any type of data where the keys are a tuple
            of tuples of of zeros and nones indicating the unique microstate.
    Returns:
        int: the number of protonation centres.

    """
    return len(tuple(data.keys())[0][0])


# DEPRECATED FUNCTIONS AND JUNK


# DEPRECATED
# def micro_constants2(molecule, lvl1, lvl2, lvl3, elvl1=None, elvl2=None,
#                      elvl3=None):
#     r"""Calculate the micro constants from parameters.
# 
#     .. math::
#         {\rm p}\bar{K}_{A\{s_j\}} = {\rm p}\bar{K}_i
#         -\sum_j{\varepsilon_{ij}s_j}
#         -frac12 \sum_{jk}{\lambda_{ijk}s_js_k}
# 
#     Parameters:
#         molecule (str): the molecule symmetry
#         lvl1 (:class:`numpy.ndarray`): a 1D array with first order parameters.
#         lvl2 (:class:`numpy.ndarray`): a 2D array with second order parameters.
#         lvl3 (:class:`numpy.ndarray`): a 3D array with third order parameters.
#         elvl1 (:class:`numpy.ndarray`, optional): a 1D array with first order
#             parameter errors.
#         elvl2 (:class:`numpy.ndarray`, optional): a 2D array with second order
#             parameter errors.
#         elvl3 (:class:`numpy.ndarray`, optional): a 3D array with third order
#             parameter errors.
# 
#     """
#     n_protcentr = len(lvl1)
#     retv = {}
#     for ms in itertools.product((0, 1), repeat=n_protcentr):
#         retv2 = {}
#         for pos, bit in enumerate(ms):
#             if bit == 1:
#                 continue
#             key = molecule[pos]
#             if key in retv2:
#                 continue
#             aux2 = sum(lvl2[pos, j]*ms[j] for j in range(n_protcentr))
#             aux3 = sum(lvl3[pos, j, k]*ms[j]*ms[k] for j in range(n_protcentr)
#                        for k in range(n_protcentr))
#             mk = lvl1[pos] - aux2 - 0.5*aux3
#             if elvl1 is None:
#                 retv2[key] = mk
#             else:
#                 eaux2 = sum(elvl2[pos, j]**2*ms[j] for j in range(n_protcentr))
#                 eaux3 = sum(elvl3[pos, j, k]**2*ms[j]*ms[k]
#                             for j in range(n_protcentr)
#                             for k in range(n_protcentr))
#                 emk = (elvl1[pos] - eaux2 - 0.5*eaux3)**0.5
#                 retv2[key] = (mk, emk)
#         retv[ms] = retv2
#     return retv


# DEPRECATED
# def microstate_multiplicity(molecule, microstate):
#     """Calculate the multiplicity of a microstate.
# 
#     Parameters:
#         molecule (str): The molecule symmetry
#         microstate (str): The string representation of a microstate
#     Returns:
#         int: the multiplicity of the microstate
#     Note:
#         The molecule and microstate strings must be expanded.
# 
#     >>> microstate_multiplicity('AABC', 'BC')
#     1
#     >>> microstate_multiplicity('AABC', 'ABC')
#     2
# 
#     """
#     count1 = collections.Counter(molecule)
#     count2 = collections.Counter(microstate)
#     remainders = (1 + count1[k] - v for k, v in count2.items())
#     return functools.reduce(lambda a, b: a*b, remainders, 1)



# DEPRECATED
# def fit_free_energy(microk, lvl1=None, lvl2=None, error=False):
#     r"""Calculate the free energy of a system based on parameters.
# 
#     The **free energy** of the system is decomposed in first, second and third
#     order parameters and calculated according to the following equation.
# 
#     .. math::
#         \frac{\beta F(\{s_j\})}{\ln 10} = -\sum_j{p\hat{K}_js_j}
#            + \frac1{2!}\sum_{ij}\varepsilon_{ij}s_is_j
#            + \frac1{3!}\sum_{ijk}\lambda_{ijk}s_is_js_k
# 
#     """
#     centres = len(microk)
#     ln10 = 2.3025851
# 
#     def zero_term(ms):
#         return sum(a*b for a, b in zip(ms, microk))
# 
#     def first_term(ms):
#         return sum(lvl1[i, j]*ms[i]*ms[j]
#                    for i, j in iprod(range(centres), range(centres)))
# 
#     def scnd_term(ms):
#         return sum(lvl2[i, j, k]*ms[i]*ms[j]*ms[k]
#                    for i, j, k in iprod(range(centres), range(centres),
#                                         range(centres)))
# 
#     def gen():
#         for microstate in iprod((0, 1), repeat=centres):
#             sum1 = zero_term(microstate)
#             if not error:
#                 sum1 *= -1
#             sum2 = first_term(microstate) if lvl1 is not None else 0.0
#             sum3 = scnd_term(microstate) if lvl2 is not None else 0.0
#             yield microstate, ln10*(sum1 + sum2/2.0 + sum3/6.0)
# 
#     return {d: k for d, k in gen()}


# DEPRECATED
# def remove_equivalent_states(molecule, data):
#     """Remove microstates that are equivalent based on molecule symmetry.
# 
#     Parameters:
#         molecule (str): a string representing the molecule symmetry
#         data (dict): any data where the keys are the microsite numerical
#             representation.
#     Returns:
#         dict: The same *data* but with the entries that are symmetrically
#             equivalente removed.
# 
#     """
#     guard = set()
#     retv = {}
#     for ms, v in data.items():
#         name = name_microstate(molecule, ms)
#         if name in guard:
#             continue
# 
#         retv[ms] = v
#         guard.add(name)
#     return retv


# def total_microstates(microstates):
#     """Calculate total number of microstates."""
#     return sum(len(m) for m in microstates)


# DEPRECATED
def name_terms(molecule, level):
    """Name the terms of the molecule at a given level.

    Parameters:
        molecule (str): A string representing the :ref:`def_molecule_symmetry`.
        level (int): the level (1, 2, or 3) of the parameters
    Returns:
        list: the sorted names of the terms
    Example:
        >>> name_terms('AABC', 1)
        ['A', 'B', 'C']
        >>> name_terms('AABC', 2)
        ['AA', 'AB', 'AC', 'BC']
        >>> name_terms('AABC', 3)
        ['AAB', 'AAC', 'ABC']

    """
    cmb = itertools.combinations(molecule, level)
    return sorted(set(("".join(q) for q in cmb)))


# DEPRECATED
# def name_microstate(molecule: str, ms: tuple) -> str:
#     """Given a numeric microstate, compute the text version of it.
# 
#     Parameters:
#         molecule (str): A string representing the :ref:`def_molecule_symmetry`.
#         ms (tuple): a 0/1 tuple representing the :ref:`microstate <microstate_definition>`.
#     Returns:
#         str: the name of the microstate.
# 
#     >>> name_microstate('AABC', (1, 0, 0, 1))
#     'AC'
#     >>> name_microstate('AAB', (1, 0, 0))
#     'A'
# 
#     """
#     return "".join(t for t, n in zip(molecule, ms) if n)


# DEPRECATED
# def molecule_info(molecule):
#     """Info given a molecule symmetry.
# 
#     Parameters:
#         molecule (str): A string representing the :ref:`def_molecule_symmetry`.
# 
#     Returns:
#         dict: info about the molecule symmetry
# 
#         * *expanded* the molecule with multipliers expanded
#         * *n_params_1* the number of first-order parameters
#         * *n_params_2* the number of second-order parameters
#         * *n_params_3* the number of third-order parameters
#         * *n_protctrs* the number of protonation centres
# 
#     """
#     info = dict()
#     expanded = expand_name(molecule)
#     info['expanded'] = expanded
#     info['n_params_1'] = len(set(expanded))
#     info['n_params_2'] = len(order_terms(expanded, 2))
#     info['n_params_3'] = len(order_terms(expanded, 3))
#     info['n_protctrs'] = len(expanded)
#     return info


# DEPRECATED
# def order_terms(molecule, n):
#     """Return valid combinations of sites of order n.
# 
#     Given a molecule, n-th order combinations of the sites defined by the
#     symmetry. Results are sorted.
# 
#     Parameters:
#         molecule (str): a string defining the summetry of the molecule.
#         n (int): the order of interactions, usually n=2 or n=3.
#     Returns:
#         list: valid n-th order combinations
# 
#     """
#     lst = ["".join(_) for _ in set(itertools.combinations(molecule, n))]
#     lst.sort()
#     return lst
