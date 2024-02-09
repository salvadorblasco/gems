# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                      !
# !                                                                             !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@protonmail.com> !
# ! Licensed under MIT license (see file LICENSE)                               !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""GUI routines for GEMS."""

import contextlib
import io
import tkinter.ttk as ttk
import tkinter as tk
import tkinter.scrolledtext
import webbrowser

import numpy as np
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NavigationToolbar

from . import __version__
import gems.libio
import gems.plotter
import gems.report
import gems.libuk


class Application(tk.Frame):
    """Main application class."""
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        master.wm_title(f'GEMS {__version__}')

        menu = self.create_menu()
        self.master.config(menu=menu)
        self.plotters = dict()
        self.create_widgets()
        self.pack()

        # TEST ONLY
        # self.text_input.delete("1.0", tk.END)
        # _file_ = '../data/py22'
        # _file_ = '../data/uvdata.dat'
        # with open(_file_, 'r') as f:
        #     self.text_input.insert("1.0", f.read())

    def create_menu(self):
        """Create all menus."""
        menu = tk.Menu(self.master)

        file_menu = tk.Menu(menu)
        file_menu.add_command(label="Open file", command=self.open_file)
        file_menu.add_command(label="Save input", command=self.save_input)
        file_menu.add_command(label="Save output", command=self.save_output)
        file_menu.add_command(label="Exit", command=self.exit_program)
        menu.add_cascade(label="File", menu=file_menu)

        options_menu = tk.Menu(menu)
        self.order = tk.IntVar(value=3)
        # options_menu.add_radiobutton(label="1st order", variable=self.order, value=1)
        # options_menu.add_radiobutton(label="1st+2nd order", variable=self.order, value=2)
        # options_menu.add_radiobutton(label="1st+3nt+3rd order", variable=self.order, value=3)

        # options_menu.add_separator()

        self.uvvis = tk.BooleanVar(value=False)
        self.dryrun = tk.BooleanVar(value=False)
        options_menu.add_radiobutton(label="NMR data", variable=self.uvvis, value=False)
        options_menu.add_radiobutton(label="UV-vis data", variable=self.uvvis, value=True)
        options_menu.add_checkbutton(label="dry run", variable=self.dryrun, onvalue=True, offvalue=False)

        menu.add_cascade(label="Options", menu=options_menu)

        help_menu = tk.Menu(menu)
        help_menu.add_command(label="Help menu", command=self._help)
        help_menu.add_command(label="About", command=self._about)
        menu.add_cascade(label="Help", menu=help_menu)

        return menu

    def create_widgets(self):
        """Create all widgets in the frame. """
        self.noteb = ttk.Notebook(self.master)
        frame_input = ttk.Frame(self.noteb)
        frame_output = ttk.Frame(self.noteb)
        self.noteb.add(frame_input, text='Input')
        self.noteb.add(frame_output, text='Output')
        self.noteb.pack(expand=1, fill=tk.BOTH)

        self.text_input = tkinter.scrolledtext.ScrolledText(frame_input, wrap=tk.NONE)
        self.text_input.pack(expand=1, fill=tk.BOTH)
        self.text_output = tkinter.scrolledtext.ScrolledText(frame_output, state='disabled',
                                                             wrap=tk.NONE)
        self.text_output.pack(expand=1, fill=tk.BOTH)

        button = tk.Button(text='Run', command=self.run_fitting)
        button.pack(side=tk.RIGHT)

    def create_navigation_toolbar(self, master, canvas):
        toolbar = NavigationToolbar(canvas, master)
        toolbar.update()
        return toolbar

    def create_plotter(self, title):
        """Create a new plotter."""
        frame = ttk.Frame(self.noteb)
        self.noteb.add(frame, text=title)
        plotter = gems.plotter.Plotter(frame)
        self.create_navigation_toolbar(frame, plotter)
        plotter.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        return plotter

    def exit_program(self):
        """Exit."""
        self.master.destroy()

    def open_file(self):
        """Open dialog to open a new file."""
        filename = tk.filedialog.askopenfilename(title="Select file")
        if not filename:
            return
        self.text_input.delete("1.0", tk.END)
        with open(filename, 'r') as fhandler:
            self.text_input.insert("1.0", fhandler.read())

    def run_fitting(self):
        """Read data and perform fitting."""
        if self.uvvis.get():
            loadf = gems.libio.load_spectra_stream
            plotf = self._plot_spectra
        else:
            loadf = gems.libio.load_stream
            plotf = self._plot

        with io.StringIO() as inputstream:
            inputstream.write(self.text_input.get('1.0', 'end'))
            inputstream.seek(0)

            try:
                data_args = loadf(inputstream)
            except:
                tk.messagebox.showerror('Error', 'Problem loading file.')
                return

            if self.uvvis.get():
                # data_args = gems.libio.load_spectra_stream(inputstream)
                title, molecule, ini_vals, keywords, pH, wavelength, yvalues = data_args
                labels = wavelength
            else:
                # data_args = gems.libio.load_stream(inputstream)
                title, molecule, ini_vals, keywords, labels, pH, yvalues = data_args

        with contextlib.redirect_stdout(io.StringIO()) as output:
            order = self.order.get()
            # infodict = gems.fit.aux_fitting1(title, pH, yvalues, molecule, ini_vals, order,
            #                                  self.uvvis.get(), False)
            infodict = gems.fit.run_fitting(title, pH, yvalues, molecule, keywords, ini_vals, order,
                                            self.uvvis.get(), self.dryrun.get())
            infodict['labels'] = labels

            if self.uvvis.get():
                infodict['wavelength'] = wavelength
                gems.fit.postfit(infodict, smooth_points=2*yvalues.shape[0])
            else:
                gems.fit.postfit(infodict)

            gems.report.print_report(infodict)
            self.text_output['state'] = tk.NORMAL
            self.text_output.delete("1.0", tk.END)
            self.text_output.insert("1.0", output.getvalue())
            self.text_output['state'] = tk.DISABLED

        plotf(infodict)

        # if self.uvvis.get():
        #     self._plot_spectra(infodict)
        # else:
        #     self._plot(infodict)

    def save_input(self):
        """Save the text in input tab to a file."""
        _save(self.text_input)

    def save_output(self):
        """Save the text in input tab to a file."""
        _save(self.text_output)

    def _about(self):
        msg = "GEMS\nThe GEneral Microspeciation Solver\n(C) by Dr Salvador Blasco 2019-2024\
               <salvador.blasco@protonmail.com>"
        tk.messagebox.showinfo(master=self, title='About GEMS', message=msg)

    def _help(self):
        """Open help in a browser."""
        docs = '../doc/build/html/index.html'
        webbrowser.open_new(docs)

    def _recall_plotter(self, label):
        """Access to a given plotter or create it if it does not exist.

        Parameters:
            label (str): the label that identifies the plotter

        The plotter is cleared in any case.
        """
        if label in self.plotters:
            plotter = self.plotters[label]
        else:
            plotter = self.create_plotter(label)
            self.plotters[label] = plotter
        plotter.clear()
        return plotter

    def _plot(self, infodict):
        labels = infodict['labels']
        group_label = {label[0] for label in labels}
        groups = dict()
        for glabel in group_label:
            plotter_data = self._recall_plotter(f'data {glabel}')
            plotter_derivs = self._recall_plotter(f'derivs {glabel}')
            lst = [n for n, lbl in enumerate(labels) if lbl[0] == glabel]
            llabels = [lbl for lbl in labels if lbl[0] == glabel]
            axes = plotter_data.figure.add_subplot(111)
            gems.plotter.plot_shifts(axes, lst, llabels, **infodict)
            axes = plotter_derivs.figure.add_subplot(111)
            gems.plotter.plot_dshifts(axes, lst, llabels, **infodict)
            groups[glabel] = lst

        self.__plot_distribution(infodict)
        self.__plot_microconstants(infodict)
        self.__plot_free_energy(infodict)
        self._update_canvases()

    def _plot_spectra(self, infodict):
        """Plot everything for spectra.

        Pops new windows with the plots for all the information contained in
        the dict.
        """
        pH = infodict['xvalues']
        absorbance = infodict['yvalues']
        wavelength = infodict['wavelength']

        plotter = self._recall_plotter('data')
        plotter.figure.clf()
        axes = plotter.figure.add_subplot(111)
        axes.plot(wavelength, absorbance.T, 'k:')
        axes.plot(wavelength, infodict['ysmooth'][::2].T, 'b-', lw=0.3)
        axes.set_xlabel('wavelength / nm')
        axes.set_ylabel('absorbance')

        plotter = self._recall_plotter('UV B-matrix')
        plotter.figure.clf()
        axes = plotter.figure.add_subplot(111)
        gems.plotter.plot_uvbmatrix(axes, **infodict)

        self.__plot_distribution(infodict)
        self.__plot_microconstants(infodict)
        self.__plot_free_energy(infodict)

        plotter = self._recall_plotter('residuals')
        plotter.figure.clf()
        axes = plotter.figure.subplots(nrows=2, ncols=1)
        y = (infodict['residuals']/infodict['yvalues'])**2
        gems.plotter.plot_uv_residuals(axes, (wavelength, pH), y)

        plotter = self._recall_plotter('weights')
        plotter.figure.clf()
        axes = plotter.figure.subplots(nrows=2, ncols=1)
        axes[0].plot(wavelength, absorbance.T, color='black')
        axes[0].set_ylabel('absorbance')
        axes[1].plot(wavelength, gems.fit.spectra_weighting(absorbance).T, color='blue')
        axes[1].set_ylabel('weight')
        axes[1].set_xlabel('wavelength')

    def __plot_distribution(self, infodict):
        plotter = self._recall_plotter('distribution')
        plotter.figure.clf()
        axes = plotter.figure.add_subplot(111)
        gems.plotter.plot_distribution(axes, **infodict)

    def __plot_free_energy(self, infodict):
        plotter = self._recall_plotter('free energy')
        plotter.figure.clf()
        axes = plotter.figure.add_subplot(111)
        gems.plotter.plot_energies(axes, **infodict)

    def __plot_microconstants(self, infodict):
        plotter = self._recall_plotter('microconstants')
        plotter.figure.clf()
        axes = plotter.figure.add_subplot(111)
        lg10_macrocts = np.log10(infodict['macroconstants'])
        stepwise_macroconstants = [a-b for a, b in zip(lg10_macrocts[1:],
                                                       lg10_macrocts[:-1])]
        gems.plotter.plot_microconstants(axes, stepwise_macroconstants, **infodict)

    def _update_canvases(self):
        for plot in self.plotters.values():
            plot.draw()


def _save(what):
    filename = tk.filedialog.asksaveasfilename(title="Select file",
                                               filetypes=(("all files", "*.*"),))
    if not filename:
        return

    data = what.get('1.0', 'end')
    with open(filename, 'w') as filehandler:
        filehandler.write(data)
