GEMS Changelog
==============

Version 0.7
-----------
* Added routines for fitting UV-vis spectra.
* Added menus to GUI for fitting options
* (lots of) Code cleanup to avoid duplicated code
* Added 'About' dialog in GUI

Version 0.6
-----------
* Code cleaning
* Bug corrected on wrong population calculations
* new :func:`gems.libuk.microstate_population`.
* Improved setup.py
* Deprecated functions removed from :mod:`gems.libuk`.
* Added :func:`gems.libio.import_numpy` to factilitate type conversion when importing from NumPy arrays.

Version 0.5
-----------

* Code cleaning
* Added option -v for choosing the verbosity of output
* Added :func:`gems.report.minimal_report` for minumum verbosity
* Added :func:`gems.libio.load_file_spectra` for loading spectrophotometric data.
* updated :func:`gems.libio.export_numpy` to used *infodict* instead of the
    deprecated *pre_plotting*.
* :func:`gems.report.report_molecule` and :func:`gems.report.generate_scheme` 
    moved to :mod:`gems.report`.
* Major internal rewriting of :mod:`gems.plotter` and :mod:`gems.fit`
* Some functions are deprecated
* Main script renaming from *gems.py* to *rungems.py*
* Submodule *gems.run* renamed *gems.report*
* Major reorganisation of *src* folder into module *gems* and script *gems.py*

Version 0.4
-----------

* Added plot :py:func:`gems.plot_energies`.
* Added Maxima Analysis
* Removed unnecesary lstsq calls for calculating B-matrix
* Added GUI with Tkinter. Created :py:mod:`gems.gui`
* moved some functions to :mod:`gems.fit` and :mod:`gems.libio`
* Many long functions split in more comprehensive chunks.
* Added help menu
* Added save input and output
* Added command line option '-w' for saving data as numpy array.

Version 0.3
-----------

* Added auxiliary function :func:`gems.libuk.num_prot_centres`
* Improved :func:`libuk.error_macro_constants`, :func:`gems.libuk.macro_constants` and :func:`gems.libuk.conditional_probability2`.


