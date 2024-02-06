# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                   !
# !                                                                          !
# ! Copyright (c) 2020-2021 by Salvador Blasco <salvador.blasco@gmail.com>   !
# ! Licensed under MIT license (see file LICENSE)                            !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""Collection of routines for plotting data.

Plotting functions with pyplot
------------------------------

* :func:`do_plot`
* :func:`plot_fitting`


Plotting functions in Axes
--------------------------

* :func:`plot_distribution`
* :func:`plot_dshifts`
* :func:`plot_microconstants`
* :func:`plot_energies`
* :func:`plot_shifts`
* :func`plot_uvbmatrix`
* :func`plot_uv_residuals`

Auxiliary plotting functions
----------------------------

* :func:`point_microconstants`
* :func:`annotate_microconstants`
* :func:`make_lines`

"""


import itertools
import math

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigureCanvas

import gems.libuk

mpl.use('TkAgg')


class Plotter(FigureCanvas):
    """Override FigureCanvas."""
    def __init__(self, master):
        self.figure = mpl.figure.Figure()
        super().__init__(self.figure, master)

    def clear(self):
        self.figure.clear()


def do_plot(infodict):
    """Plot everything for NMR data.

    Pops new windows with the plots for all the information contained in
    the dict.
    """

    plot_fitting(**infodict)

    plot_distribution(pl.gca(), **infodict)
    pl.figure()

    lg10_macrocts = np.log10(infodict['macroconstants'])
    stepwise_macroconstants = [a-b for a, b in zip(lg10_macrocts[1:],
                                                   lg10_macrocts[:-1])]
    plot_microconstants(pl.gca(), stepwise_macroconstants, **infodict)
    pl.figure()

    plot_energies(pl.gca(), **infodict)
    pl.show()


def do_plot_spectra(infodict):
    """Plot everything for spectra.

    Pops new windows with the plots for all the information contained in
    the dict.
    """
    pH = infodict['xvalues']
    absorbance = infodict['yvalues']
    wavelength = infodict['wavelength']

    pl.plot(wavelength, absorbance.T, 'k:')
    pl.plot(wavelength, infodict['ysmooth'][::2].T, 'b-', lw=0.3)
    pl.xlabel('wavelength / nm')
    pl.ylabel('absorbance')
    pl.figure()

    # plot_fitting(**infodict)
    plot_uvbmatrix(pl.gca(), **infodict)
    pl.figure()

    plot_distribution(pl.gca(), **infodict)
    pl.figure()

    lg10_macrocts = np.log10(infodict['macroconstants'])
    stepwise_macroconstants = [a-b for a, b in zip(lg10_macrocts[1:], lg10_macrocts[:-1])]
    plot_microconstants(pl.gca(), stepwise_macroconstants, **infodict)
    pl.figure()

    plot_energies(pl.gca(), **infodict)

    fig = pl.figure()
    axs = fig.subplots(nrows=2, ncols=1)
    y = (infodict['residuals']/infodict['yvalues'])**2
    gems.plotter.plot_uv_residuals(axs, (wavelength, pH), y)

    fig, axs = pl.subplots(nrows=2, ncols=1, sharex=True)
    axs[0].plot(wavelength, absorbance.T, color='black')
    axs[0].set_ylabel('absorbance')
    axs[1].plot(wavelength, gems.fit.spectra_weighting(absorbance).T, color='blue')
    axs[1].set_ylabel('weight')
    axs[1].set_xlabel('wavelength')

    pl.show()


def plot_distribution(axes, xsmooth, distribution, centres_occupation,
                      centres_unique, n_protctrs, **kwargs):
    """Plot the macrospecies and microspecies distribution.

    Parameters:
        axes (:class:`matplotlib.Axes.axes`): the axes to plot into.
        xsmooth:
        distribution:
        centres_occupation:
        centres_unique
        n_protctrs
        kwargs (dict): Ignored

    .. seealso:: :func:`gems.fit.postfit`
    """
    axes.plot(xsmooth, distribution, '-')
    axes.plot(xsmooth, centres_occupation, ':', linewidth=2.5)
    axes.legend([str(a) for a in range(n_protctrs+1)] + centres_unique)
    axes.set_xlabel('pD')
    axes.set_ylabel('population')


def plot_uv_residuals(axes, xdim, residuals):
    """Plot residuals.

    Parameters:
        axes (:class:`matplotlib.Axes.axes`): the axes to plot into.
    """
    for n, (x_, txt) in enumerate(zip(xdim, ('wavelength', 'pH'))):
        axes[n].stem(x_, np.sum(residuals, axis=n), use_line_collection=True)
        axes[n].set_xlabel(txt)
        axes[n].set_ylabel('residual')
        axes[n].set_yscale('log')


def plot_uvbmatrix(axes, wavelength, bmatrix, **kwargs):
    ydata = bmatrix.T
    axes.plot(wavelength, ydata)
    axes.legend([str(i) for i in range(ydata.shape[1])])
    axes.set_xlabel('wavelength / nm')
    axes.set_ylabel('absorbance')


def plot_shifts(axes, lst, llabels, xvalues, yvalues, xsmooth, ysmooth,
                distribution, centres_occupation, **kwargs):
    """Plot chemical shifts and distribution in background.
    """
    axes.plot(xsmooth, ysmooth[:, lst], color='grey')
    for lbl, y in zip(llabels, yvalues[:, lst].T):
        axes.plot(xvalues, y, 'o', label=lbl)
    axes.legend()
    axes.set_xlabel('pD')
    axes.set_ylabel('$\\delta$ / ppm')
    ax2 = axes.twinx()
    ax2.set_ylabel('population')
    ax2.plot(xsmooth, distribution, '-', linewidth=0.5)
    ax2.plot(xsmooth, centres_occupation, ':')
    ax2.set_navigate(False)


def plot_dshifts(axes, lst, llabels, xsmooth, ysmooth, centres_occupation, centres_unique,
                 **kwargs):
    axes.plot(xsmooth, -np.gradient(ysmooth[:, lst], xsmooth, axis=0))
    axes.legend(llabels, loc='lower right')
    axes.axhline(lw=0.5, color='black')
    axes.set_xlabel('pD')
    axes.set_ylabel(r'$-{\rm d}\delta / {\rm dpD}$')
    twin_ax = axes.twinx()
    twin_ax.plot(xsmooth, -np.gradient(centres_occupation, xsmooth, axis=0), ':')
    twin_ax.legend(centres_unique, loc='lower left')
    twin_ax.set_ylabel(r'$-{\rm d}\theta/{\rm dpD}$')
    twin_ax.set_navigate(False)


def plot_energies(axes, molecule, free_energy, **kwargs):
    """Plot energies.
    """
    for ms, energy in free_energy.items():
        n = sum(ms)
        name = gems.libuk.name_microstate(molecule, ms)
        axes.annotate(name, (n, energy), xytext=(0, 10),
                      textcoords='offset points', ha='center')
        axes.hlines(energy, n-0.4, n+0.4, colors='blue', label='X')

    axes.set_ylabel(r'$\Delta G/RT$')
    axes.set_xlabel('degree of protonation')
    axes.set_xticks(range(len(molecule)+1))


def plot_microconstants(axes, stepwise_macroconstants, msnames, microconstants,
                        n_protctrs, **kwargs):
    """Plot microconstants."""
    axes.scatter(np.arange(1, n_protctrs + 1), stepwise_macroconstants,
                 s=50.0, color='red', zorder=3.0)

    out_x, out_y, labels = point_microconstants(msnames, microconstants)

    axes.scatter(out_x, out_y, color='blue', s=30.0, zorder=2.0)
    for line in make_lines(microconstants):
        l2d = mpl.lines.Line2D(line[0], line[1], linewidth=0.2,
                               color='grey', zorder=1.5)
        axes.add_line(l2d)
    annotate_microconstants(out_x, out_y, labels, axes)
    statk = [math.log10(k) + stepwise_macroconstants[0]
             for k in gems.libuk.kdeflation(n_protctrs)]
    axes.plot(np.arange(1, n_protctrs + 1), statk, color='black',
              zorder=1.0, linestyle=':')
    axes.set_xlabel('Degree of protonation, $n$')
    axes.set_ylabel(r'$\log_{10} K$')
    axes.set_xticks(np.arange(1, n_protctrs + 1))


def plot_fitting(xvalues, yvalues, xsmooth, ysmooth, distribution,
                 centres_occupation, labels, centres_unique, **kwargs):
    """Plot fitting."""

    group_label = {label[0] for label in labels}
    groups = dict()
    for glabel in group_label:
        lst = [n for n, lbl in enumerate(labels) if lbl[0] == glabel]
        llabels = [lbl for lbl in labels if lbl[0] == glabel]
        plot_shifts(pl.gca(), lst, llabels, xvalues, yvalues, xsmooth, ysmooth,
                    distribution, centres_occupation)
        pl.figure()
        plot_dshifts(pl.gca(), lst, llabels, xsmooth,
                     ysmooth, centres_occupation,
                     centres_unique)
        pl.figure()
        groups[glabel] = lst


# TODO rename this function to make it internal to the module
def make_lines(micro_constants):
    """Calculate where to draw the energy levels.

    Parameters:
        micro_constants (dict): the microconstants
    Returns:
        list: the coordinates of the lines to be drawn in the format
            ((x1, y1), (x2, y2))

    """
    keys = micro_constants.keys()
    lines = []
    for ua, ub in itertools.product(keys, keys):
        a = ua[0]
        b = ub[0]
        if sum(a) + 1 == sum(b):
            diff = (abs(i-j) for i, j in zip(b, a))
            if sum(diff) == 1:
                mk1 = micro_constants[ua]
                mk2 = micro_constants[ub]
                common_keys = set(mk1.keys()) & set(mk2.keys())
                for k in common_keys:
                    lines.append(((1+sum(a), 1+sum(b)), (mk1[k], mk2[k])))
    return lines


# TODO rename this function to make it internal to the module
def annotate_microconstants(xvalues, yvalues, labels, axes):
    """Insert microconstants annotation into graph.

    Parameters:
        xvalues (iterable): the x values of the label position.
        yvalues (iterable): the y values of the label position.
        labels (iterable): the labels
        axes (:class:`matplotlib.axes.Axes`): the axes to plot into.

    """
    properties = {'xytext': (10.0, 10.0), 'textcoords': 'offset points',
                  'arrowprops': {'arrowstyle': 'wedge'},
                  'va': 'bottom', 'ha': 'left'}

    for x, y, l in zip(xvalues, yvalues, labels):
        axes.annotate(l, (x, y), **properties)


def point_microconstants(msnames, micro_constants):
    """Calculate positions for microconstant labels.

    Parameters:
        molecule (str):
        micro_constants (iterable):
    Returns:
        tuple: where [0] and [1] are the *x*  and *y*  coordinates of the label,
            and [2] are the labels.

    """
    out_x = []
    out_y = []
    labels = []
    # breakpoint()
    for ms, msteps in micro_constants.items():
        _root_label_ = msnames[ms]
        # _root_label_ = gems.libuk.name_microstate(molecule, ms)
        if _root_label_ == '':
            _root_label_ = 'Ø'
        n = sum(ms[0]) + 1
        for mstep, mk in msteps.items():
            added = msnames[mstep].replace(_root_label_, '')
            new_label = _root_label_ + '+' + added
            if new_label in labels:
                continue

            out_x.append(n)
            out_y.append(mk)
            labels.append(new_label)
    return out_x, out_y, labels
