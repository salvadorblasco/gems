# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""This module contains the reading/writing and import/export routines.

Import/export fitting data to/from NumPy array
* :func:`import_numpy`
* :func:`export_numpy`

Load data from file
* :func:`load_stream`
* :func:`load_spectra_stream`
* :func:`load_file_spectra`
* :func:`load_file`

Auxiliary/Internal functions
* :func:`test_modules`
* :func:`print_welcome`
* :func:`prune_data`
"""

import re
import sys
import numpy as np
import scipy


def export_numpy(filename: str, infodict: dict):
    """Write all information into NumPy format file.

    Parameters:
        filename (:class:`str`): the file name.
        infodict (:class:`dict`): all the information (see :func:`gems.fit.postfit`)

    """
    np.savez_compressed(filename, **infodict)


def import_numpy(filename: str):
    """Read information from NumPy format file.

    Uses :func:`numpy.load` and then applies some type conversions.

    Parameters:
        filename (str): the file to read from.
    Returns:
        dict: the data read

    .. versionadded:: 0.6
    """
    with np.load(filename, allow_pickle=True) as fhandler:
        _i0 = lambda a: a.item(0)
        _fmap = {'molecule': str, 'result': _i0, 'order': int, 'labels': tuple, 'message': str,
                 'free_energy': _i0, 'error_free_energy': _i0, 'microstate_probability': _i0,
                 'conditional_probability': _i0, 'microconstants': _i0, 'centres_unique': list}
        result = {k: fhandler[k] if k not in _fmap else _fmap[k](fhandler[k]) for k in fhandler}
    return result


def load_file(filename: str):
    """Open file and load its contents.

    A text file containing the data to be fit. See
    :ref:`data_file_format`.

    Arguments:
        filename (:class:`str`): The file to read from
    Returns:
        tuple: title (str), molecule (str), initial values (list of floats),
            labels (list of strings), pH (:class:`numpy.ndarray`),
            shifts (:class:`numpy.ndarray`).

    .. seealso:: :func:`gems.libio.loadstream`
    """
    with open(filename, 'r') as fhandler:
        retv = load_stream(fhandler)

    return retv


def load_file_spectra(filename: str):
    """Open file containing spectral data and load its contents.

    Arguments:
        filename (str): The file to read from
    Returns:
        tuple: title (str), molecule (str), initial values (list of floats),
            pH (:class:`numpy.ndarray`),
            wavelength (:class:`numpy.ndarray`),
            absorbances (:class:`numpy.ndarray`).

    .. seealso:: :func:`load_file`
    """
    with open(filename, 'r') as fhandler:
        retv = load_spectra_stream(fhandler)

    return retv


def load_spectra_stream(streamio):
    """Read datastream containing spectral data.

    Arguments:
        streamio (StreamIOr): The data stream
    Returns:
        tuple: title (str), molecule (str), initial values (list of floats),
            pH (:class:`numpy.ndarray`),
            wavelength (:class:`numpy.ndarray`),
            absorbances (:class:`numpy.ndarray`).

    .. seealso:: :func:`load_file_spectra`
    """
    header = __load_header(streamio)
    pH = np.array([float(s) for s in streamio.readline().split()])
    d = np.genfromtxt(streamio, usemask=True, missing_values='X', filling_values=0.0)
    wavelength = d[:, 0]
    absorbances = d[:, 1:].T

    return *header, pH, wavelength, absorbances


def load_stream(streamio):
    """Load data from data stream.

    For details on the file data format, see :ref:`data_file_format`.

    Arguments:
        streamio (file): The strem to read data from.
    Returns:
        tuple: title (str), molecule (str), initial values (list of floats),
            labels (list of strings), pH (:class:`numpy.ndarray`),
            shifts (:class:`numpy.ndarray`).

    """
    header = __load_header(streamio)
    labels = streamio.readline().split()

    d = np.genfromtxt(streamio, missing_values='X', usemask=True, filling_values=0.0)
    np.ma.masked_invalid(d)

    _ = d[:, 0]
    if np.ma.is_masked(_):
        raise RuntimeError("Missing values not allowed in pH/pD column")

    pH = np.copy(_.data)

    _ = d[:, 1:]
    if np.ma.is_masked(_):
        shifts = np.ma.copy(_)
        shifts.fill_value = 0.0
    else:
        shifts = np.copy(_)

    return *header, labels, pH, shifts


def prune_data(array, rows_ignored, cols_ignored, row_labels, col_labels):
    """Remove unused data from the bulk of data.

    Parameters:
        array (:class:`numpy.ndarray`): the initial data.
        rows_ignored (sequence):
        cols_ignored (sequence):
        row_labels (sequence):
        col_labels (sequence):

    Returns:
        tuple: item[0] is :class:`numpy.ndarray` the pruned data, item[1] and item[2] are
            the respective row_labels and col_labels pruned.
    """
    pruned_data = np.delete(np.delete(array, rows_ignored, axis=0),
                            cols_ignored, axis=1)
    row_pruned = np.delete(row_labels, rows_ignored)
    col_pruned = np.delete(col_labels, cols_ignored)
    return pruned_data, row_pruned, col_pruned


def print_welcome():
    """Print welcome message."""
    welcome = """
                        GGGG   EEEEE  MM   MM   SSSS
                       G       E      M M M M  S
                       G  GGG  EEE    M  M  M   SSS
                       G    G  E      M     M      S
                        GGGG   EEEEE  M     M  SSSS

                    The GEneral Microspeciation Solver
                   (C) by Dr Salvador Blasco 2019-2024
                     <salvador.blasco@protonmail.com>


    Discaimer: This is experimental software, provided "as is". No warranty
    is provided. You are not excused from critically analyse the output.

    """
    print(welcome)


# def process_keywords(args):
#     rekw = re.compile(r'(fix|restrain) (.*)', re.IGNORECASE)
#     retv = []
#     for arg in args:
#         if (m := rekw.fullmatch(arg.strip())) is None:
#             print("Warning: keyword not recognised. It will be ignored")
#             print(" > ", arg)
#         else:
#             retv.append((m.group(1), m.group(2)))
#     return retv


def test_modules():
    """Test the modules."""
    print('running:')
    print('  * python  {0}.{1}.{2}'.format(*sys.version_info))
    print('  * numpy  ', np.version.version)
    print('  * scipy  ', scipy.version.version)
    print()


def __load_header(streamio):
    title = streamio.readline().strip()
    keywords = {}
    while True:
        line = streamio.readline().strip().upper()
        match line.split():
            case ['$MATRIX', var1]:
                n = int(var1)
                mt = " ".join(streamio.readline().strip() for _ in range(n))
                m = np.fromstring(mt, dtype=int, sep=' ').reshape((n,n))
                keywords['matrix'] = m
            case ['$FIX', param]:
                __push_param(keywords, 'fix', param)
            case ['$RESTRAIN', *param]:
                __push_param(keywords, 'restrain', tuple(param))
            case ['$BOUND', variable, lower_bound, upper_bound]:
                __push_param(keywords, 'bound', (variable, float(lower_bound), float(upper_bound)))
            case other:
                if line[0] == '$':
                    raise ValueError(f"Invalid keyword: {_}")
                break
    molecule = line
    aux = streamio.readline().strip()
    ini_vals = [float(s) for s in aux.split()]
    return title, molecule, ini_vals, keywords


def __push_param(mydict: dict, kw: str, param):
    if kw in mydict:
        mydict[kw].append(param)
    else:
        mydict[kw] = [param]
