# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                   !
# !                                                                          !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@gmail.com>   !
# ! Licensed under MIT license (see file LICENSE)                            !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""Routines for printing information from the data fitting.

Public API

* :func:`gems.report.report_molecule` for printing a report based on the molecule symmetry.
* :func:`gems.report.minimal_report` for a small report
* :func:`gems.report.print_report` for a full report. It includes the minimal report.
* :func:`gems.report.generate_scheme`

Auxiliary functions

* :func:`gems.report.section`
* :func:`gems.report.print_residuals`
* :func:`gems.report.print_nuclei_influence`
* :func:`gems.report.print_microstates`
* :func:`gems.report.print_deltamatrix`
* :func:`gems.report.print_correlation`
* :func:`gems.report.print_bmatrix`
* :func:`gems.report.print_amatrix`
* :func:`gems.report.print_array_half`
* :func:`gems.report.print_array_full`
* :func:`gems.report.print_array_sorted`
* :func:`gems.report.print_parameters`
* :func:`gems.report.print_populations`
* :func:`gems.report.prefit`

"""

import enum
import itertools

import numpy as np

import gems.fit
import gems.libuk


class Verbosity(enum.Enum):
    VERBOSITY_MINIMAL = enum.auto()
    VERBOSITY_REGULAR = enum.auto()
    VERBOSITY_VERBOSE = enum.auto()


verbosity = Verbosity.VERBOSITY_REGULAR


def prefit(title, n_points, n_data):
    """Print initial message."""
    print('file successfully loaded')
    print('TITLE: ', title)
    print(f'  {n_points} titration points')
    print(f'  {n_data} data points')
    print()


def minimal_report(infodict):
    """Print minimal report."""
    print()
    if infodict['success']:
        print('Optimization successful. {message}'.format(**infodict))
    else:
        print('Optimization failure. {message}'.format(**infodict))

    print('iterations: ', infodict['iterations'])
    print('evaluations: ', infodict['function_evaluations'])
    print('final chi squared : ', infodict['chisq'])
    print()

    print_parameters(infodict)

    if infodict['correlation'] is not None:
        print_correlation(**infodict)
    print()

    if 'wavelength' in infodict:
        print_residuals_uv(infodict['residuals'], infodict['xvalues'],
                           infodict['yvalues'], infodict['wavelength'])
    else:
        print_residuals(infodict['residuals'], infodict['xvalues'],
                        infodict['yvalues'], infodict['labels'])


def print_parameters(infodict: dict):
    nms = infodict['msnames']
    parms = infodict['parameters']
    mapping = infodict['mapping']
    section('parameters')
    ordinals = ('first', 'second', 'third')
    for n, ordn in enumerate(ordinals):
        print(f"{ordn:>12}-order:")
        for ums in gems.libuk.filter_by_macro(mapping, 1+n):
            name = nms[ums]
            ms0 = gems.libuk.pop_microstate(ums)
            val = parms.vmc(ms0)
            err = parms.evmc(ms0)
            stat = parms.microstate_status(ms0)
            print(f'      {name:>8}:  {val:>10.4f} ± {err:.4f}  {stat}')

def print_report(infodict):
    """Print report.

    Parameters:
        infodict (dict):  a :class:`dict` containing all the infomation.

    .. seealso: :func:`fit.postfit`

    """
    print_molecule(**infodict)
    if infodict['result'] is not None:
        minimal_report(infodict)
    else:
        print_parameters(infodict)

    print_microstates(infodict)
    print_macroconstants(**infodict)
    print_bmatrix(**infodict)
    print_amatrix(**infodict)
    print_deltamatrix(**infodict)
    if 'wavelength' not in infodict:
        print_nuclei_influence(**infodict)
    print_populations(**infodict)


def print_populations(n_microcts, n_protctrs, centres_unique, idx_unique,
                      theta, xvalues, macrostate_probability, **kwargs):
    """Print information of population analysis.

    Parameters:
        n_microcts (int): the number of microconstants in the system
        n_protctrs (int): the number of protonation centres
        centres_unique
        idx_unique
        theta
        xvalues
        macrostate_probability
        kwargs (dict): this is ignored

    """
    _title1 = 'Degree of occupation'.center(n_microcts * 8, ' ')
    _title2 = 'Population'.center((n_protctrs+1) * 8, ' ')
    print(10*' ' + _title1 + 5*' ' + _title2)
    print('    pD' + 10*' ' + (7*' ').join(centres_unique) + '   ' +
          ((n_protctrs + 1)*'{:9}').format(*range(n_protctrs+1)))
    fmt = '{:8.2f}' + len(idx_unique) * ' {:8.3f}' + ' |' + \
          (n_protctrs+1) * ' {:8.3f}'
    for ph_, th, cmapr in zip(xvalues,
                              n_protctrs * theta[idx_unique, :].T,
                              macrostate_probability):
        print(fmt.format(ph_, *th, *cmapr))


# def print_parameters(fit1, err1, fit2, err2, fit3, err3, molecule):
#     """Print report about the parameters."""
#     def _faux(ttl, fit, err, lbls):
#         txt = '      {:>8}:  {:>10.4f} ± {:.4f}'
#         print('  {}-order:'.format(ttl))
#         for f, e, t in zip(fit, err, lbls):
#             print(txt.format(t, f, e))
# 
#     section('refined parameters')
#     _faux('first', fit1, err1, gems.libuk.name_terms(molecule, 1))
#     if fit2 is not None:
#         _faux('second', fit2, err2, gems.libuk.name_terms(molecule, 2))
#     if fit3 is not None:
#         _faux('third', fit3, err3, gems.libuk.name_terms(molecule, 3))


def print_array_sorted(array, labels, n=10):
    """Print the *n* highest values of the array.

    Parameters:
        array (:class:`numpy.ndarray`): The array to print
        labels (sequence): the labels for each row/column of the array
        n (int, optional): the number of entries to print

    """
    import itertools
    lng = len(array)
    d = {(labels[i], labels[j]): array[i, j]
         for i, j in itertools.product(range(lng), range(lng)) if i < j}
    print(' term1    term2       r')
    print(' -------+--------+----------')
    sorted_data = sorted(d, key=lambda x: abs(d[x]), reverse=True)
    for k in itertools.islice(sorted_data, n):
        fmt = '  {:<8} {:<8} {:8.4f}'
        print(fmt.format(k[0], k[1], d[k]))


def print_array_full(array, labels):
    """Print a full array with labels.

    Parameters:
        array (:class:`numpy.ndarray`): The array to print
        labels (sequence): the labels for each row/column of the array

    """
    lng = len(array)
    fmt1 = 12*' ' + '{:^10}' * lng
    print(fmt1.format(*labels))
    fmt2 = '{:10}' + lng*'{:10.4f}'
    for i in range(len(array)):
        print(fmt2.format(labels[i], *array[i]))


def print_array_half(array, labels):
    """Print a triangular array with labels.

    Parameters:
        array (:class:`numpy.ndarray`): The array to print
        labels (sequence): the labels for each row/column of the array

    """
    lng = len(array)
    for i in range(1, lng):
        print('{:10}'.format(labels[i]), end='')
        for j in range(lng):
            if i > j:
                print('{:10.4f}'.format(array[i, j]), end='')
            else:
                print(5*' ', end='')
        print()
    fmt1 = 12*' ' + '{:^10}' * (lng-1)
    print(fmt1.format(*labels))


def print_amatrix(amatrix, molecule, **kwargs):
    """Print report on :ref:`def_a_matrix`.

    Parameters:
        amatrix (:class:`numpy.ndarray`): the A matrix.
        molecule (str): A string representing the expanded molecule symmetry.
        kwargs (dict): this is ignored

    """
    n_protctrs = amatrix.shape[1] - 1
    fmt = '{:8}' + (n_protctrs + 1) * ' {:10.3f}'
    section('A matrix')
    _header = (n_protctrs + 1) * '{:^11d}'
    print('      n = ', _header.format(*range(n_protctrs + 1)))
    guard = set()
    for site, aline in zip(molecule, amatrix):
        if site in guard:
            continue
        else:
            print(fmt.format(site, *aline))
            guard.add(site)


def print_bmatrix(bmatrix, labels, **kwargs):
    """Print report on :ref:`def_b_matrix`.

    Parameters:
        bmatrix
        labels
        kwargs (dict): this is ignored

    """
    n_protctrs = len(bmatrix) - 1
    fmt = '{:8}' + (n_protctrs + 1) * ' {:10.3f}'
    section('B matrix')
    _header = (n_protctrs + 1) * '{:^11d}'
    print('      n = ', _header.format(*range(n_protctrs + 1)))
    for lbl, b in zip(labels, bmatrix.T):
        print(fmt.format(lbl, *b))


def print_correlation(correlation, parameters, msnames, **kwargs):
    """Print information on the correlation matrix."""
    section('correlation')
    sorted_params = parameters.get_sorted_params()
    name_params = [msnames[i] for i in sorted_params]
    if len(name_params) > 6:
        print_array_sorted(correlation, name_params)
    else:
        print_array_half(correlation, name_params)


def print_deltamatrix(dmatrix, delta_span, labels, molecule, **kwargs):
    """Print report on :ref:`def_delta_matrix`.

    Parameters:
        dmatrix
        delta_span
        labels
        molecule
        kwargs (dict): this is ignored

    """
    centres_unique = sorted(set(molecule))
    idx_unique = [molecule.index(site) for site in centres_unique]
    print()
    print(' Delta matrix '.center(30, '-').center(60))
    fmt = '{:10}' + len(idx_unique) * ' {:10.3f}' + ' - {:10.3f}'
    view_delta_matrix = dmatrix[:, idx_unique]
    print(17*' ' + (10*' ').join(centres_unique) + 11*' ' + 'span')
    for lbl, d, s in zip(labels, view_delta_matrix, delta_span.T):
        print(fmt.format(lbl, *d, s))
    print()


def print_macroconstants(log_macroconstants, error_log_macroconstants, **kwargs):
    """Print the macroconstant report.

    Parameters:
        log_macroconstants (:class:`numpy.ndarray`):
        error_log_macroconstants (:class:`numpy.ndarray`):
        kwargs (dict): this is ignored

    """
    section('Calculated macroconstants')
    print('  n      log10K')
    print(' ---   ----------')
    lbeta = log_macroconstants[1:]
    errlbeta = error_log_macroconstants[1:]
    for n, (b, e) in enumerate(zip(lbeta, errlbeta), start=1):
        print('  {}     {:>.3f} ± {:<.3f}'.format(n, b, e))
    print()
    print("    If these constants represent deuteration equilibria, the analogous")
    print(" protonation equilibrium can be estimated with the relation:")
    print("              pKa(D2O) = 0.32 + 1.044 pKa(H2O)")
    print("              [R. Delgado, et al. Anal. Chim. Acta 1991, 245, 271-282]")


def print_microstates(infodict: dict):
    """Print the microstates.

    Parameters:
        molecule (str): A string representing the expanded molecule symmetry.
        free_energy (dict): the free energy
        error_free_energy (dict): the error of the free energy
        conditional_probability (dict): the conditional microstate probability
        microconstants (iterable):
        kwargs (dict): this is ignored

    """
    molecule: str = infodict['molecule']
    msteps: dict = infodict['microsteps']
    nms: dict = infodict['msnames']
    mapping: dict = infodict['mapping']
    condprob: dict = infodict['conditional_probability']
    free_energy: dict = infodict['free_energy']
    error_free_energy: dict = infodict['error_free_energy']
    # ums_population = gems.libuk.merge_equivalent_microstates(mapping, condprob, lambda x, y: x+y)
    # ums_free = gems.libuk.merge_equivalent_microstates(mapping, free_energy, lambda x, y: x)
    # ums_efree = gems.libuk.merge_equivalent_microstates(mapping, error_free_energy, lambda x, y: x)
    ums_population = infodict['ums_population']
    ums_free = infodict['ums_free_energy']
    ums_efree = infodict['ums_error_free_energy']
    microconstants = infodict['microconstants']

    section('Microstates and Microsteps Analysis')
    for level in range(1+len(molecule)):
        print(f'n = {level}')
        for ums in gems.libuk.filter_by_macro(mapping, level):
            name = nms[ums]
            if not name:
                name = '∅'
            popu = ums_population[ums]
            free = ums_free[ums]
            efree = ums_efree[ums]
            print(f'{name:>5s} -> F = {free:6.2f} ± {efree:6.2f}  pop = {100*popu:6.2f}')
            if level < len(molecule):
                for mstp in msteps[ums]:
                    # prod = nms[mstp]
                    ms0 = gems.libuk.pop_microstate(ums)
                    mstp0 = gems.libuk.pop_microstate(mstp)
                    site = "".join(_m for _m, _a, _b in  zip(molecule, mstp0, ms0) if _a != _b)
                    #muk = gems.libuk.micro_constant(ums, mstp, ums_free)
                    muk = microconstants[ums][mstp]
                    print(f'      + {site:10} {muk:10.3f}')

    # sorted_keys = sorted(free_energy.keys(), key=sum)
    # level = -1
    # done = set()
    # for key in sorted_keys:
    #     if sum(key) > level:
    #         level = sum(key)
    #         print('n = ', level)
    #     name = gems.libuk.name_microstate(molecule, key)
    #     if name in done:
    #         continue

    #     done.add(name)
    #     name_ = name if name != '' else 'Ø'
    #     err = error_free_energy[key]**0.5
    #     nwe = '{:.2f} ± {:.2f}'
    #     fmt = '    {:5} -> F = {:15} pop = {:<10.2f}'
    #     fnwe = nwe.format(free_energy[key], err)
    #     print(fmt.format(name_, fnwe, 100*microstate_population[name]))
    #     for mstep, mk in microconstants[key].items():
    #         print('      + {:10} {:10.3f}'.format(mstep, mk))


def print_molecule(input_molecule, molecule, connectivity, isomorphisms, keywords={}, **kwargs):
    """Print report of a given molecule.

    Parameters:
        molecule (str): A string representing the expanded molecule symmetry.

    """
    print(f"Input molecule: {input_molecule}")
    print()
    print("Expanded molecule: ", molecule, '\n')
    if 'matrix' in keywords:
        print("Given connectivity matrix")
        repr_connectivity(connectivity, molecule)
    else:
        print("Tentative connectivity matrix")
        repr_connectivity(connectivity, molecule)
        print("   NOTE: if this matrix does not correctly represent your molecule you have to input the matrix by hand.")
    print()
    print(f"found {len(isomorphisms)} isomorphisms")
    print()


def print_nuclei_influence(xsmooth, ysmooth, labels, centres_unique, centres_occupation,
                           **kwargs):
    """Print the influence of each nucleus in each protonation step.

    Parameters:
        xsmooth (:class:`numpy.ndarray`):
        ysmooth (:class:`numpy.ndarray`):
        labels (sequence):
        centres_unique (:class:`numpy.ndarray`):
        centers_occupation (:class:`numpy.ndarray`):
        kwargs (dict): this is ignored

    """
    section('Maxima analysis')
    import scipy.signal
    gc = np.gradient(centres_occupation, axis=0)
    gs = np.gradient(ysmooth, axis=0)
    maxs = scipy.signal.argrelmax(gc**2)

    for idx, col in np.nditer(maxs):
        print(f'pD = {xsmooth[idx]:.2f}',
              'protonation' if gc[idx, col] < 0 else 'deprotonation',
              f'on {centres_unique[col]}')
        dct = {lbl: gs[idx, n] for n, lbl in enumerate(labels)}

        print('  Nucleus    -1000dδ/dpD')
        for k in sorted(dct, key=lambda k: (k[0], abs(dct[k])), reverse=True):
            print(f'  {k:8}  {-1000*dct[k]:>8.2f}')
        print()
    print()


def print_residuals(residuals, pD, shifts, labels, num=20):
    """Print the part of the report dealing with residuals.

    Given a list of residuals produced by the fitting, print the larges
    **n** values.

    Parameters:
        residuals (:class:`numpy.ndarray`):
        pD (:class:`numpy.ndarray`): the x-values
        shifts (:class:`numpy.ndarray`): the y-values
        labels (iterable): the labels associated with y-values
        num (int, default=20): the number of residuals to print

    """
    # _num_ = min(num, shifts.size)
    # section('worst points')
    # shape = residuals.shape
    # indices = np.squeeze(
    #     np.dstack(
    #         np.unravel_index(
    #             np.argsort(np.abs(residuals.ravel())), shape)))[::-1]
    # print('  #  point    pD      shift  residual   nucleus')
    # fmt = '{:3d} {:4d} {:8.2f} {:9.3f} {:9.3f}   {}'

    # count = 0
    # for r, c in indices:
    #     if shifts[r, c] is np.ma.masked:
    #         continue
    #     print(fmt.format(count + 1, 1 + r, pD[r], shifts[r, c], residuals[r, c],
    #                      labels[c]))
    #     count += 1
    #     if count > _num_:
    #         break

    header = '  #  point    pD      shift  residual   nucleus'
    _print_residuals(residuals, pD, shifts, labels, header, num)


def print_residuals_uv(residuals, pH, absorbance, wavelength, num=20):
    """Print the part of the report dealing with residuals.

    Given a list of residuals produced by the fitting, print the larges
    **n** values.

    Parameters:
        residuals (:class:`numpy.ndarray`):
        pD (:class:`numpy.ndarray`): the x-values
        shifts (:class:`numpy.ndarray`): the y-values
        labels (iterable): the labels associated with y-values
        num (int, default=20): the number of residuals to print

    """
    # _num_ = min(num, absorbance.size)
    # section('worst points')
    # shape = residuals.shape
    # indices = np.squeeze(
    #     np.dstack(
    #         np.unravel_index(
    #             np.argsort(np.abs(residuals.ravel())), shape)))[::-1]
    # print('  #  point   pH       absorb.  residual wavelength')
    # fmt = '{:3d} {:4d} {:8.2f} {:9.3f} {:9.3f}   {}'

    # count = 0
    # for r, c in indices:
    #     if absorbance[r, c] is np.ma.masked:
    #         continue
    #     print(fmt.format(count + 1, 1 + r, pH[r], absorbance[r, c], residuals[r, c], wavelength[c]))
    #     count += 1
    #     if count > _num_:
    #         break

    header = '  #  point   pH       absorb.  residual wavelength'
    _print_residuals(residuals, pH, absorbance, wavelength, header, num)



def _print_residuals(residuals, xvalues, yvalues, lvalues, header, num=20):
    """Print the part of the report dealing with residuals.

    Given a list of residuals produced by the fitting, print the larges
    **n** values.

    Parameters:
        residuals (:class:`numpy.ndarray`):
        pD (:class:`numpy.ndarray`): the x-values
        shifts (:class:`numpy.ndarray`): the y-values
        labels (iterable): the labels associated with y-values
        num (int, default=20): the number of residuals to print

    """
    _num_ = min(num, yvalues.size)
    section('worst points')
    shape = yvalues.shape
    indices = np.squeeze(
        np.dstack(
            np.unravel_index(
                np.argsort(np.abs(residuals.ravel())), shape)))[::-1]
    print(header)
    fmt = '{:3d} {:4d} {:8.2f} {:9.3f} {:9.3f}   {}'

    count = 0
    for r, c in indices:
        if yvalues[r, c] is np.ma.masked:
            continue
        print(fmt.format(count + 1, 1 + r, xvalues[r], yvalues[r, c], residuals[r, c], lvalues[c]))
        count += 1
        if count > _num_:
            break


def generate_scheme(molecule):
    """Generate the possible protonation pathways for a given molecule.

    Parameters:
        molecule (str): A string representing the expanded molecule symmetry.

    Returns:
        list: each element is a :class:`dict` whose keys are the possible
            microstates and the values are sets which contain the possible
            protonation sites available.

    >>> gems.report.generate_scheme('AAB')                                           
    [{'': {'A', 'B'}},
     {'B': {'A'}, 'A': {'A', 'B'}},
     {'AB': {'A'}, 'AA': {'B'}},
     {'AAB': set()}]

    """
    n_protctrs = len(molecule)

    output = [{} for _ in range(n_protctrs+1)]
    for microstate in itertools.product((False, True), repeat=n_protctrs):
        n = sum(microstate)
        di = output[n]
        name = gems.libuk.name_microstate(molecule, microstate)
        if name not in di:
            microsteps = {molecule[pos]
                          for pos, bit in enumerate(microstate)
                          if not bit}
            di[name] = microsteps
    return output


def report_molecule(input_molecule):
    """Print report of a given molecule.

    Parameters:
        molecule (str): A string representing the expanded molecule symmetry.

    """
    print(f"Input molecule: {input_molecule}")
    print()
    molecule = gems.libuk.expand_name(input_molecule)
    print("Expanded molecule: ", molecule, '\n')
    connect = gems.libuk.compute_connectivity(molecule)
    print("Tentative connectivity matrix")
    repr_connectivity(connect, molecule)
    print("   NOTE: if this matrix does not correctly represent your molecule you have to input the matrix by hand.")
    print()
    isomorphisms = gems.libuk.compute_isomorphisms(molecule, connect)
    print(f"found {len(isomorphisms)} isomorphisms")
    print()

    mapping = gems.libuk.compute_mapping(molecule, isomorphisms)
    names = gems.libuk.name_equivalent_microstates(molecule, mapping)
    total_microstates = 0
    parameters = [0, 0, 0]
    for level in range(1+len(molecule)):
        print(5*'-', 'macrostate n =', level, 5*'-')
        print("m: n:           equivalent microstates")
        for ums in gems.libuk.filter_by_macro(mapping, level):
            total_microstates += 1
            if 0 < level < 4:
                parameters[level-1] += 1
            print(f"{len(ums):2d} {names[ums]:<12}", ", ".join(ms_to_text(molecule, ms) for ms in ums)) 
        print()

    print(f"Total of {total_microstates} independent microstates")
    print()
    print("Independent parameters to fit:", sum(parameters))
    for tx, nn in zip((' first', 'second', ' third'), parameters):
        print(f"  {tx}-level: {nn}")

    # molec = gems.libuk.expand_name(molecule)
    # scheme = generate_scheme(molec)
    # print()
    # total_microstates = 0
    # total_microsteps = 0
    # print('molecule symmetry: ', molecule, '\n')
    # print('molecule expanded: ', molec, '\n')
    # print()
    # for n, microlevel in enumerate(scheme):
    #     n_mstates = len(microlevel)
    #     total_microstates += n_mstates
    #     txt_mists = ", ".join(gems.libuk.compact_name(_)
    #                           for _ in microlevel.keys())
    #     if txt_mists == '':
    #         txt_mists = 'Ø'
    #     print(f'n = {n}, {n_mstates} microstate(s): ', txt_mists)
    #     n_microsteps = 0
    #     for microstate, microsteps in microlevel.items():
    #         if microstate == '':
    #             _microstate = 'Ø'
    #         else:
    #             _microstate = gems.libuk.compact_name(microstate)
    #         for microstep in microsteps:
    #             aux = "".join(sorted(microstate + microstep))
    #             _microresult = gems.libuk.compact_name(aux)
    #             print('  ', _microstate, ' + ', microstep, ' -> ',
    #                   _microresult)
    #         n_microsteps += len(microsteps)
    #     if n < len(scheme) - 1:
    #         print(f'  {n_microsteps} microstep(s)')
    #     print()
    #     total_microsteps += n_microsteps
    # print(f'total microstates = {total_microstates}')
    # print(f'total microsteps = {total_microsteps}')

    # print()
    # print('Fitting information')
    # n_microl = len(set(molec))
    # n_term2 = len(gems.libuk.order_terms(molec, 2))
    # n_term3 = len(gems.libuk.order_terms(molec, 3))
    # print('  first order parameters:  ', n_microl)
    # print('  second order parameters: ', n_term2)
    # print('  third order parameters:  ', n_term3)
    # print('  - - - - - - - - - - - - - -')
    # total_params = n_microl + n_term2 + n_term3
    # print('  total parameters:        ', total_params)


def section(text):
    """Separate sections."""
    w = max(30, len(text) + 10)
    print()
    print((' ' + text + ' ').center(w, '-').center(60))
    print()


def repr_connectivity(connectivity, molecule):
    print("  ", "".join(f"{m:2}" for m in molecule))
    for rlet, row in zip(molecule, connectivity):
        print(f"{rlet}", "".join(f"{x:2}" for x in row))

def repr_clasified_microstates(my_set):
    print('{', end='')
    for t in my_set:
        print('(', end='')
        print(",".join(('(' + "".join(str(_) for _ in tt) + ')' for tt in t)), end='')
        print(')', end=',')
    print('}')


def ms_to_text(molecule, ms):
    return "".join(molecule[n].upper() if b else molecule[n].lower() for n, b in enumerate(ms))


def print_isomorphisms(isomorphisms, molecule):
    for n, isomorphism in enumerate(isomorphisms):
        print(f'isomorphism #{n}')
        repr_connectivity(isomorphism, molecule)
