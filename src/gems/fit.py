# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""All the routines related to fitting.

Public Functions:
* :func:`do_fit`
* :func:`calc_error_free_energy`
* :func:`fit_shifts2`
* :func:`calc_residuals2`
* :func:`fit_shifts`
* :func:`calc_residuals`
* :func:`rearrange_parameters`
* :func:`arrange_parameters`

Auxiliary Functions:
* :func:`postfit`
* :func:`fobj`
* :func:`first_term`
* :func:`build_terms`
"""


import collections
import sys

import numpy as np
import scipy.optimize

import gems.libuk
import gems.vmc


def run_fitting(title, pH, yvalues, molecule, keywords, ini_vals, order, flag_uv, dry_run=False):
    gems.report.prefit(title, len(pH), yvalues.size)
    expmolec = gems.libuk.expand_name(molecule)

    if 'matrix' in keywords:
        connectivity = keywords['matrix']
    else:
        connectivity = gems.libuk.compute_connectivity(expmolec)
    isomorphisms = gems.libuk.compute_isomorphisms(expmolec, connectivity)
    mapping = gems.libuk.compute_mapping(expmolec, isomorphisms)
    msteps = gems.libuk.microsteps(len(expmolec), mapping)
    ms_naming = gems.libuk.name_equivalent_microstates(expmolec, mapping)
    parameters = gems.vmc.FittingParams(mapping, ms_naming)
    parameters.set_values(ini_vals)
    _apply_keywords(keywords, parameters, ms_naming)

    proton_activity = 10**-pH
    if flag_uv:
        weights = spectra_weighting(yvalues)
    else:
        weights = 1.0

    print('Starting...')

    # method = 'L-BFGS-B'
    # method = 'BFGS'
    # method='Nelder-Mead'
    method = keywords.get('method', 'BFGS')

    print(f' optimize with {method}')
    # print(parameters.get_bounds())
    # for ums in parameters.get_sorted_params():
    #     print(ums, ms_naming[ums])
    # print(parameters.parameters)
    # print(parameters.ids)
    # sys.exit(1)
    # _bounds = ((None, None), (None, None), (None, None), (2.0, 4.0), (None, None))
    _bounds = parameters.get_bounds()

    if _bounds is not None and method != 'Nelder-Mead':
        print(" WARNING: bounds not allowed for method {method}. Bounds will be ignored.")
        _bounds = None

    fobj = fobj_preparation(parameters, order, proton_activity, yvalues, weights, calc_residuals)
    chisq0 = fobj(ini_vals)
    print(f' initial chi squared: {chisq0}\n')

    if dry_run:
        retv = None
    else:
        retv = scipy.optimize.minimize(fobj, ini_vals, method=method, bounds=_bounds)
    infodict = {'xvalues': pH, 'yvalues': yvalues, 'result': retv,
                'input_molecule': molecule,
                'molecule': expmolec, 'order': order, 'parameters': parameters,
                'mapping': mapping, 'isomorphisms': isomorphisms, 'connectivity': connectivity,
                'microsteps': msteps, 'msnames': ms_naming}
    return infodict


def fobj_preparation(parameters,
                     order: int, 
                     proton_activity: np.ndarray,
                     yvalues: np.ndarray,
                     weights: np.ndarray,
                     fresiduals: callable):
    """Generate the objective function.
    """
    def fobj(values):
        # print(values)
        parameters.set_values(values)
        free_enrg = _calc_free_energy(parameters)
        residual = fresiduals(free_enrg, yvalues, proton_activity)
        chisq = np.sum(weights*residual**2)
        return chisq
    return fobj


def _apply_keywords(args, parameters, ms_names):
    for kw in args:
        match kw:
            case 'fix':
                fixes = args[kw]
                for fix in fixes:
                    parameters.create_constraint(fix)
            case 'restrain':
                restraints = args[kw]
                for restraint in restraints:
                    parameters.create_restraint(restraint)
            case 'matrix':
                pass
            case 'method':
                pass
            case 'bound':
                _process_bound(args[kw], parameters, ms_names)
            case other:
                raise ValueError(f"Invalid keyword: {other}")


def _process_bound(bounds, parameters, ms_naming):
    reverse_names = {v: k for k, v in ms_naming.items()}
    for umsn, lower, upper in bounds:
        ums = reverse_names[umsn]
        parameters.create_bound(ums, lower, upper)

def _calc_free_energy(parameters):
    microstates = gems.libuk.generate_all_microstates(parameters.size)
    return {ms: parameters.free_energy(ms) for ms in microstates}


def _calc_free_energy_error(parameters):
    microstates = gems.libuk.generate_all_microstates(parameters.size)
    return {ms: parameters.error_free_energy(ms) for ms in microstates}


# def arrange_parameters(expmolec, initvals, order=3):
#     """Build the initial parameter array from initial guess.
# 
#     Construct an initial array filled with zeros and then insert
#     the given parameters in its proper places, leaving the not given
#     parameters as zeros. Additional initial values given are ignored.
# 
#     Parameters:
#         expmolec (str): a string defining the expanded symmetry of the
#             molecule.
#         initvals (iterable): list of initial values
#         order (int): the order of interactions (1, 2 or 3)
#     Returns:
#         :class:`numpy.ndarray`: A 1D array with the parameters arranged in
#             the proper order.
#     Raises:
#         ValueError: if the order parameter is not within 1--3.
# 
#     """
#     if not 0 < order <= 3:
#         raise ValueError
#     n_params_1 = len(set(expmolec))
#     n_params_2 = len(gems.libuk.order_terms(expmolec, 2))
#     n_params_3 = len(gems.libuk.order_terms(expmolec, 3))
#     n = (n_params_1, n_params_2, n_params_3)
#     total_params = sum(n[:order])
# 
#     parameters = np.zeros(total_params)
#     for i in range(min(len(parameters), len(initvals))):
#         parameters[i] = initvals[i]
#     return parameters


# def rearrange_parameters(param, param_errors, molecule, order):
#     """Rearrange a flattened array of parameters into order-sorted arrays.
# 
#     Parameters:
#         param (iterable): a 1D array containing ordered parameters
#         param_errors (iterable): a 1D array containing ordered parameter errors
#         molecule (str): A string representing the molecule symmetry
#         order (int): A number indicating the interaction order level (1..3)
# 
#     Returns:
#         tuple: fit_params1, err_params1, fit_params2, err_params2,
#             fit_params3, err_params3
# 
#     """
#     n_params_1 = len(set(molecule))
#     n_params_2 = len(gems.libuk.order_terms(molecule, 2))
#     _lim = n_params_1 + n_params_2
#     fit_params1 = param[:n_params_1]
#     err_params1 = param_errors[:n_params_1]
# 
#     fit_params2, err_params2, fit_params3, err_params3 = 4*(None,)
# 
#     if order >= 2:
#         fit_params2 = param[n_params_1:_lim]
#         err_params2 = param_errors[n_params_1:_lim]
# 
#     if order == 3:
#         fit_params3 = param[_lim:]
#         err_params3 = param_errors[_lim:]
# 
#     return fit_params1, err_params1, fit_params2, err_params2, fit_params3, \
#         err_params3


# DEPRECATED
# def do_fit(init_vars, order, expmolec, ah, shifts, weights=1.0):
#     """Fit data.
# 
#     Calls :func:`scipy.optimization.minimize` in order to minimize
#     :func:`fobj`.
# 
#     Parameters:
#         init_vars (:class:`numpy.ndarray`):
#         order (int):
# 
#     .. versionchanged:: 0.5
#         Added option *weights*
# 
#     """
#     method = 'BFGS'
#     # method='Nelder-Mead'
#     chisq0 = fobj(init_vars, order, expmolec, ah, shifts, weights)
#     res = scipy.optimize.minimize(fobj, init_vars, method=method,
#                                   args=(order, expmolec, ah, shifts, weights))
#     return chisq0, res


# DEPRECATED
# def build_terms(molecule, args=None):
#     """Construct epsilon and lambda matrices from parameters.
# 
#     Parameters:
#         molecule (str): a string defining the summetry of the molecule.
#         args (sequence): a sequence containing the parameters
# 
#     Returns:
#         tuple: the first element is a 2D array containing the binary terms;
#             the second is a 3D array containing the tertiary terms.
# 
#     """
#     if args is None or len(args) == 0:
#         return None, None
# 
#     def two_equal(seq):
#         count = collections.Counter(seq)
#         return any(i > 1 for i in count.values())
# 
#     neps = len(molecule)
# 
#     array2 = np.zeros((neps, neps))
#     terms2 = gems.libuk.order_terms(molecule, 2)
#     params2 = args[:len(terms2)]
#     dump = [(2, array2, terms2, params2)]
# 
#     array3 = None
#     if len(args) > len(terms2):
#         array3 = np.zeros((neps, neps, neps))
#         terms3 = gems.libuk.order_terms(molecule, 3)
#         params3 = args[len(terms2):]
#         dump.append((3, array3, terms3, params3))
# 
#     for n, array, terms, params in dump:
#         it = np.nditer(array, op_flags=['writeonly'], flags=['multi_index'])
#         while not it.finished:
#             txt = "".join(sorted(molecule[_] for _ in it.multi_index))
#             if two_equal(it.multi_index):
#                 it.iternext()
#                 continue
#             term = terms.index(txt)
#             it[0] = params[term]
#             it.iternext()
# 
#     return array2, array3


# def first_term(molecule, terms):
#     """Return first term."""
#     run = sorted(set(molecule))
#     return [terms[run.index(site)] for site in molecule]


# DEPRECATED
# def fobj(x, *args):
#     """Objective function to be used with :func:`minimize`.
#
#     Parameters:
#         x (:class:`numpy.ndarray`): A 1D array with the fitting parameters
#         args[0] (int): the parameter range to fit (1, 2 or 3).
#         args[1] (str): the molecule symmetry.
#         args[2] (:class:`numpy.ndarray`): the activity of protons
#         args[3] (:class:`numpy.ndarray`): the experimental data where the
#             titration points are arranged in rows and nuclei are arranged
#             in columns.
#         args[4] (float or :class:`numpy.ndarray`): the weights to be used for
#             minimization.
#     Returns:
#         float: the sum of the residuals squared.
#
#     """
#     fit_range, molecule, proton_activity, shifts, weights = args
#     n_microcts = len(set(molecule))
#     # n_protctrs = len(molecule)
#     lgpk = first_term(molecule, x[:n_microcts])
#     term2, term3 = None, None
#     if fit_range > 1:
#         term2, term3 = build_terms(molecule, x[n_microcts:])
#     free_enrg = gems.libuk.fit_free_energy(lgpk, term2, term3)
#     _fresiduals = calc_residuals
#     # _fresiduals = calc_residuals2
#     residual = _fresiduals(free_enrg, shifts, proton_activity)
#     chisq = np.sum(weights*residual**2)
#     return chisq


def calc_residuals(free_energy, shifts, proton_activity):
    r"""Calculate the residuals.

    The fitting is done by equation
    :math:`\delta_l = \delta_l^{(0)} + \sum_{n=1}^N B_{ln} P_n (a_H)`
    by fitting the experimental chemical shifts to the calculated
    :math:`P_n(a_H)`

    Parameters:
        free_energy (:class:`dict`): the free energy
        shifts (:class:`numpy.ndarray`): the chemical shifts
        proton_activity (:class:`numpy.ndarray`): the proton activity

    Returns:
        :class:`numpy.ndarray`: the residuals

    """
    calc_beta = np.array(gems.libuk.macro_constants(free_energy))
    macro_prob = gems.libuk.macrostate_probability(calc_beta, proton_activity)
    # !!
    # yvalues = shifts
    # size = len(calc_beta) - 1
    # micro_prob = gems.libuk.microstate_probability(free_energy, proton_activity)
    # conditional_probability = gems.libuk.conditional_probability2(free_energy)
    # theta = gems.libuk.avg_prtdeg(micro_prob)
    # delta = yvalues.T @ np.linalg.pinv(theta)
    # matrix_a = gems.libuk.matrix_a(conditional_probability, size)
    # matrix_b = delta @ matrix_a
    # ycalc = macro_prob @ matrix_b.T
    # residual = ycalc - yvalues
    # !!

    ## bmatrix = fit_shifts(macro_prob, shifts)
    ## residual = np.dot(macro_prob, bmatrix) - shifts
    # calc_shifts = fit_shifts(macro_prob, shifts)
    matrix_b = gems.libuk.matrix_b(macro_prob, shifts)
    # calc_shifts = macro_prob @ matrix_b
    calc_shifts = np.dot(macro_prob, matrix_b)
    residual = calc_shifts - shifts
    # print(np.sum(residual**2))
    return residual


# DEPRECATED. USe gems.libuk.matrix_b
# def fit_shifts(macro_prob, shifts):
#     r"""Return calculated chemical shifts.
#
#     The calculated chemical shifts are computed according to the equation
#     below. The parameter *shifts* can be a :class:`numpy.ma.ndarray`; in this
#     case, the missing values are zeroed.
#
#     .. math::
#         \delta^{\rm calc} = P^{+}_n(a_{\mathrm{H}}) \cdot
#         \delta^{\rm exp} \cdot P_n(a_{\mathrm{H}})
#
#     Parameters:
#         macro_prob (:class:`numpy.ndarray`): The  macrostate probability.
#         shifts (:class:`numpy.ndarray`): The experimental chemical shifts.
#
#     Return:
#         :class:`numpy.ndarray`: the calculated chemical shifts.
#
#     """
#     improb = np.linalg.pinv(macro_prob)
#     if np.ma.is_masked(shifts):
#         # calc_shifts = np.dot(macro_prob, np.ma.dot(improb, shifts, strict=False))
#         # bmatrix = np.ma.dot(improb, shifts).data
#         # bmatrix = gems.libuk.mlstsq(macro_prob, shifts).data
#         # calc_shifts = np.dot(macro_prob, bmatrix)
#         bmatrix = np.ma.dot(improb, shifts)
#         calc_shifts = np.ma.dot(macro_prob, bmatrix)
#     else:
#         # bmatrix = np.dot(improb, shifts)
#         # calc_shifts = np.dot(macro_prob, np.dot(improb, shifts))
#         bmatrix = improb @ shifts
#         calc_shifts = macro_prob @ bmatrix
#     return calc_shifts
#     # return bmatrix
#     # import pudb
#     # pudb.set_trace()
#     # return np.dot(macro_prob, np.ma.dot(improb, shifts, strict=False))


def calc_residuals2(free_energy, shifts, proton_activity):
    r"""Calculate the residuals.

    The fitting is done by equation
    :math:`\delta_l = \delta_l^{(0)} + \sum_{m=1}^N \Delta_{lm} P_n (a_H)`
    by fitting the experimental chemical shifts to the calculated
    :math:`P_n(a_H)`

    Parameters:
        free_energy (:class:`dict`): the free energy
        shifts (:class:`numpy.ndarray`): the chemical shifts
        proton_activity (:class:`numpy.ndarray`): the proton activity

    Returns:
        :class:`numpy.ndarray`: the residuals

    """
    # calc_beta = np.array(gems.libuk.macro_constants(free_energy))
    # macro_prob = gems.libuk.macrostate_probability(calc_beta, ah)
    # n = len(tuple(free_energy)[0])
    msprob = gems.libuk.microstate_probability(free_energy, proton_activity)
    theta = gems.libuk.avg_prtdeg(msprob)
    calc_shifts = fit_shifts2(theta, shifts)
    residual = calc_shifts - shifts
    return residual


def fit_shifts2(theta, shifts):
    r"""Return calculated chemical shifts.

    The calculated chemical shifts are computed according to the equation
    below. The parameter *shifts* can be a :class:`numpy.ma.ndarray`; in this
    case, the missing values are zeroed.

    .. math::
        \delta^{\rm calc} = \theta \cdot
        \left(\delta^{\rm exp} \cdot \theta^+\right)

    Parameters:
        macro_prob (:class:`numpy.ndarray`): The  macrostate probability.
        shifts (:class:`numpy.ndarray`): The experimental chemical shifts.

    Return:
        :class:`numpy.ndarray`: the calculated chemical shifts.

    """
    imtheta = np.linalg.pinv(theta)
    return np.dot(theta.T, np.ma.dot(imtheta.T, shifts, strict=False))


# DEPRECATE
def calc_error_free_energy(molecule, variance, n_microcts):
    """Calculate the error in free energy.

    Parameters:
        molecule (str): the symmetry of the molecule
        variance (sequence): the variance in the free energy
        n_microcts (int): the number of microconstants

    """
    term1 = first_term(molecule, variance[:n_microcts])
    term2, term3 = build_terms(molecule, variance[n_microcts:])
    return gems.libuk.fit_free_energy(term1, term2, term3, error=True)


def spectra_weighting(data, weight_param=0.01):
    """Calculate automatic spectra weighting.

    Parameters:
        data (:class:`numpy.ndarray`):
        w0 (float):
    Returns:
        :class:`numpy.ndarray`:
    """
    return weight_param/(weight_param + np.gradient(data, axis=1)**2)


def postfit(retval, smooth_points=100):
    """Complete calculations after the fitting has been done.

    The input is a :py:class:`dict` with some data (see below). This function
    updates the dict with a lot of additional data.

    * n_nuclei (:class:`int`): the number of nuclei
    * molecule (str): the symmetry of the molecule
    * n_microcts (int): the number of microconstants in the system
    * n_protctrs (int): the number of protonation centres
    * ah (:class:`numpy.ndarray`):
    * iterations (int): the number of iterations
    * function_evaluations (int): the number of times the :func:`fobj` has
        been evaluated.
    * success
    * message
    * chisq (float): the final sum of squared residuals
    * covariance
    * fit_params1 (:class:`numpy.ndarray`):
    * fit_params2 (:class:`numpy.ndarray`):
    * fit_params3 (:class:`numpy.ndarray`):
    * err_params1 (:class:`numpy.ndarray`):
    * err_params2 (:class:`numpy.ndarray`):
    * err_params3 (:class:`numpy.ndarray`):
    * correlation
    * free_energy (dict):
    * error_free_energy
    * residuals
    * microstate_probability
    * theta
    * conditional_probability
    * macroconstants
    * error_macroconstants
    * macrostate_probability
    * microstate_population (dict):
    * bmatrix (:class:`numpy.ndarray`):
    * microconstants
    * centres_unique
    * idx_unique
    * amatrix
    * dmatrix
    * delta_span
    * smooth_ah
    * xsmooth
    * ysmooth
    * distribution
    * centres_occupation

    Parameters:
        retval (dict): a :class:`dict` that must contain the following fields
            *xvalues* (:class:`numpy.ndarray`),
            *yvalues* (:class:`numpy.ndarray`),
            *result* (:class:`scipy.optimize.OptimizationResult`),
            *labels* (sequence), *molecule* (str). This variable is modified
            during the call.
        smooth_points (int): the number of points that will be used in the
            smooth fitting result.

        * *n_nuclei*: The number of nuclei

    """
    result = retval['result']
    molecule = retval['molecule']

    n_nuclei = retval['yvalues'].shape[1]
    retval['n_nuclei'] = n_nuclei

    retval['n_microcts'] = len(set(molecule))
    retval['n_protctrs'] = len(molecule)

    retval['ah'] = 10**-retval['xvalues']
    if result is not None:
        retval['iterations'] = result['nit']
        retval['function_evaluations'] = result['nfev']
        retval['success'] = result['success']
        retval['message'] = result['message']
        retval['chisq'] = result['fun']

    if retval['result'] is not None:
        _aux_postfit1(retval)
    _aux_postfit2(retval)
    _aux_postfit3(retval)
    _aux_postfit4(retval)
    _aux_postfit5(retval, smooth_points)


def _aux_postfit1(retval):
    result = retval['result']
    # order = retval['order']
    x = result['x']
    total_params = len(x)

    if 'hess_inv' in result:
        if isinstance(result['hess_inv'], scipy.optimize._lbfgsb_py.LbfgsInvHessProduct):
            hess_inv = result['hess_inv'].todense()
        else:
            hess_inv = result['hess_inv']
        covariance = hess_inv*result['fun']/total_params
    else:
        covariance = np.zeros((total_params, total_params))
    retval['covariance'] = covariance

    param_variance = np.diag(covariance)
    param_errors = np.sqrt(param_variance)

    # params = gems.fit.rearrange_parameters(x, param_errors, retval['molecule'],
    #                                        order)
    # retval['fit_params1'], retval['err_params1'], retval['fit_params2'], \
    # retval['err_params2'], retval['fit_params3'], retval['err_params3'] = params

    if 'hess_inv' in result:
        correlation = covariance/(param_errors[:, None]*param_errors[None, :])
    else:
        correlation = None
    retval['correlation'] = correlation

    parameters = retval['parameters']
    parameters.set_values(x)
    covariance = retval['covariance']
    param_variance = np.diag(covariance)
    parameters.set_errors(np.sqrt(param_variance))


def _aux_postfit2(retval):
    parms = retval['parameters']
    # mapping = retval['mapping']

    # term1 = gems.fit.first_term(molecule, fit_params1)
    # term2, term3 = gems.fit.build_terms(molecule, x[n_microcts:])
    # free_enrg = gems.libuk.fit_free_energy(term1, term2, term3)
    # error_free_energy = gems.fit.calc_error_free_energy(molecule,
    #                                                     param_variance,
    #                                                     n_microcts)
    # ufree_energy = {ms: parms.free_energy(ms[0]) for ms in mapping.values()}
    # uerror_free_energy = {ms: parms.error_free_energy(ms[0]) for ms in mapping.values()}

    free_energy = _calc_free_energy(parms)
    error_free_energy = _calc_free_energy_error(parms)
    retval['free_energy'] = free_energy
    retval['error_free_energy'] = error_free_energy
    # retval['ufree_energy'] = ufree_energy
    # retval['uerror_free_energy'] = uerror_free_energy


def _aux_postfit3(retval):
    free_enrg = retval['free_energy']
    proton_activity = retval['ah']

    residuals = gems.fit.calc_residuals(free_enrg, retval['yvalues'], proton_activity)
    mistprob = gems.libuk.microstate_probability(free_enrg, proton_activity)
    theta = gems.libuk.avg_prtdeg(mistprob)
    condprob = gems.libuk.conditional_probability2(free_enrg)
    # mistpop = gems.libuk.microstate_population(condprob, retval['molecule'])
    calc_beta = np.array(gems.libuk.macro_constants(free_enrg))
    calc_beta_error = \
        np.array(gems.libuk.error_macro_constants(free_enrg,
                                                  retval['error_free_energy']))
    calc_macro_prob = gems.libuk.macrostate_probability(calc_beta, proton_activity)


    mapping = retval['mapping']
    free_energy = retval['free_energy']
    error_free_energy = retval['error_free_energy']
    ums_population = gems.libuk.merge_equivalent_microstates(mapping, condprob, lambda x, y: x+y)
    ums_free = gems.libuk.merge_equivalent_microstates(mapping,free_energy, lambda x, y: x)
    ums_efree = gems.libuk.merge_equivalent_microstates(mapping, error_free_energy, lambda x, y: x)

    retval['residuals'] = residuals
    retval['microstate_probability'] = mistprob
    # retval['microstate_population'] = mistpop
    retval['theta'] = theta
    retval['conditional_probability'] = condprob
    retval['macroconstants'] = calc_beta
    retval['error_macroconstants'] = calc_beta_error
    retval['macrostate_probability'] = calc_macro_prob
    retval['ums_free_energy'] = ums_free
    retval['ums_error_free_energy'] = ums_efree
    retval['ums_population'] = ums_population

def _aux_postfit4(retval):
    n_protctrs = retval['n_protctrs']
    molecule = retval['molecule']
    calc_beta = retval['macroconstants']
    calc_beta_error = retval['error_macroconstants']
    condprob = retval['conditional_probability']

    bmatrix = gems.libuk.matrix_b(retval['macrostate_probability'],
                                  retval['yvalues'])
    assert bmatrix.shape == (retval['n_protctrs'] + 1, retval['n_nuclei'])

    micro_constants = gems.libuk.micro_constants(retval['microsteps'],
                                                 retval['ums_free_energy'])
    retval['bmatrix'] = bmatrix
    retval['microconstants'] = micro_constants
    err_lg_beta = 0.434294*np.sqrt(calc_beta_error)/calc_beta
    retval['log_macroconstants'] = np.log10(calc_beta)
    retval['error_log_macroconstants'] = err_lg_beta

    centres_unique = sorted(set(molecule))
    retval['centres_unique'] = centres_unique
    idx_unique = [molecule.index(site) for site in centres_unique]
    retval['idx_unique'] = idx_unique

    delta0 = bmatrix[0, :]
    matrix_a = gems.libuk.matrix_a(condprob, n_protctrs)
    matrix_delta = np.dot(bmatrix.T - delta0[:, None],
                          np.linalg.pinv(matrix_a))
    retval['amatrix'] = matrix_a
    retval['dmatrix'] = matrix_delta

    delta_span = np.max(bmatrix, axis=0) - np.min(bmatrix, axis=0)
    retval['delta_span'] = delta_span


def _aux_postfit5(retval, smooth_points):
    xvalues = retval['xvalues']
    calc_beta = retval['macroconstants']
    free_enrg = retval['free_energy']
    n_protctrs = retval['n_protctrs']
    bmatrix = retval['bmatrix']
    idx_unique = retval['idx_unique']

    xsmooth = np.linspace(np.min(xvalues), np.max(xvalues), smooth_points)
    ah2 = 10**-xsmooth
    distribution = gems.libuk.macrostate_probability(calc_beta, ah2)
    mistprob = gems.libuk.microstate_probability(free_enrg, ah2)
    smooth_theta = gems.libuk.avg_prtdeg(mistprob)
    ysmooth = np.dot(distribution, bmatrix)
    centres_occupation = n_protctrs*smooth_theta[idx_unique, :].T

    retval['smooth_ah'] = ah2
    retval['xsmooth'] = xsmooth
    retval['ysmooth'] = ysmooth
    retval['distribution'] = distribution
    retval['centres_occupation'] = centres_occupation


# DEPRECATED use run_fitting
def aux_fitting1(title, pH, yvalues, molecule, ini_vals, order, flag_uv):
    """Auxiliary function for fitting routines."""
    gems.report.prefit(title, len(pH), yvalues.size)
    expmolec = gems.libuk.expand_name(molecule)
    assert isinstance(pH, np.ndarray)
    proton_activity = 10**-pH
    if flag_uv:
        weights = spectra_weighting(yvalues)
    else:
        weights = 1.0
    initial_values = arrange_parameters(expmolec, ini_vals, order)
    print('Starting...')
    chisq0, retv = do_fit(initial_values, order, expmolec, proton_activity, yvalues, weights)
    print(f'initial chi squared: {chisq0}')
    infodict = {'xvalues': pH, 'yvalues': yvalues, 'result': retv,
                'molecule': expmolec, 'order': order}
    return infodict
