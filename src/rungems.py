#!/usr/bin/env python3

# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                   !
# !                                                                          !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@gmail.com>   !
# ! Licensed under MIT license (see file LICENSE)                            !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""Main script of GEMS software.

It may be executed as a GUI, without arguments of in a terminal.
"""


import argparse
import sys

import gems.libuk
import gems.report
import gems.libio
import gems.plotter


__version__ = '0.9'


def _create_parser():
    str_descr = f'GEMS: The GEneral Microspeciation Solver v{__version__}'
    parser = argparse.ArgumentParser(description=str_descr)
    parser.add_argument('-s', '--symmetry',
                        help="display information on the molecule",
                        metavar='MOLECULE', type=str)
    parser.add_argument('-f', '--fit',
                        help="load data in file and start computation",
                        metavar='FILE', type=str)
    parser.add_argument('-d', '--dry-run',
                        help="load data in file and compute but do not fit",
                        default=False, action='store_true')
    parser.add_argument('-p', '--plot', help="display window with plots",
                        default=False, action='store_true')
    parser.add_argument('-o', '--order',
                        help="order of interactions to fit", metavar='ORDER',
                        type=int, choices=[1, 2, 3], default=3)
    parser.add_argument('-u', '--uvvis',
                        help="the input file contains UV-vis data",
                        default=False, action='store_true')
    parser.add_argument('-v', '--verbosity',
                        help="change output verbosity", metavar='V',
                        type=int, choices=[1, 2], default=2)
    parser.add_argument('-w', '--write',
                        help="write results in numpy array for later use",
                        metavar='FILE', type=str)
    return parser


def run():
    """Initialize program execution."""

    # Parse arguments
    parser = _create_parser()
    parsed_args = parser.parse_args()

    gems.libio.print_welcome()
    print('version: ', __version__)
    gems.libio.test_modules()

    if parsed_args.symmetry:
        gems.report.report_molecule(parsed_args.symmetry)

    if parsed_args.fit:
        order = parsed_args.order
        if parsed_args.uvvis:
            floader = gems.libio.load_file_spectra
        else:
            floader = gems.libio.load_file

        try:
            file_contents = floader(parsed_args.fit)
        except ValueError as verr:
            print("ERROR: input file contains errors")
            print(f"-> {verr}")
            print("program terminated")
            sys.exit(1)

        if parsed_args.uvvis:
            # file_contents = gems.libio.load_file_spectra(parsed_args.fit)
            title, molecule, ini_vals, keywords, pH, wavelength, yvalues = file_contents
            labels = wavelength
        else:
            # file_contents = gems.libio.load_file(parsed_args.fit)
            title, molecule, ini_vals, keywords, labels, pH, yvalues = file_contents

        infodict = gems.fit.run_fitting(title, pH, yvalues, molecule, keywords, ini_vals, order,
                                        parsed_args.uvvis, parsed_args.dry_run)
        infodict['labels'] = labels

        if parsed_args.uvvis:
            infodict['wavelength'] = wavelength
            gems.fit.postfit(infodict, smooth_points=2*yvalues.shape[0])
        else:
            gems.fit.postfit(infodict)

        infodict['keywords'] = keywords
        if parsed_args.verbosity == 2:
            gems.report.print_report(infodict)
        else:
            gems.report.minimal_report(infodict)

        if parsed_args.write:
            gems.libio.export_numpy(parsed_args.write, infodict)

        if parsed_args.plot:
            if parsed_args.uvvis:
                gems.plotter.do_plot_spectra(infodict)
            else:
                gems.plotter.do_plot(infodict)

    print('\nProgram terminated successfully')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        import tkinter as tk
        import gems.gui
        root = tk.Tk()
        root.option_add('*tearOff', False)          # remove dashed line on top of menus
        app = gems.gui.Application(master=root)
        app.mainloop()
    else:
        run()
