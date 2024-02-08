.. gems documentation master file, created by
   sphinx-quickstart on Fri Nov  8 16:34:19 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


##################
GEMS Documentation
##################

:program:`GEMS` is the **GE**\neral  **M**\icrospeciation **S**\olver. It is aimed at solving
acid-base microspeciation equilibria from NMR and spectroscopic data. It is

* Open Source
* Multiplatform
* Free
* Fast

Installing and Running GEMS
===========================

Installing
**********

:program:`GEMS` is a python script which needs numpy and scipy to run. Optionally you can
have matplotlib for graphs.

You will need

* `python 3 <http://www.python.org/>`_
* `numpy <http://www.numpy.org/>`_
* `scipy <http://www.scipy.org/>`_
* Optionally `matplotlib <http://www.matplotlib.org/>`_


First run
*********

If you execute the command :command:`gems.py -h` the usage text will show

.. literalinclude:: console_options


Fitting NMR data
****************

When run with the **-f FILE** option, the contents of the file **FILE** is loaded
and fed into the program and fitted. The contents should conform to the format 
specified in :ref:`data_file_format`.


.. _data_file_format:

Data File Format
================

GEMS can load text files containing data in the following format:

* LINE1: A title. This line is an identifier. It is ignored.
* LINE2: The symmetry of the molecule.
* LINE3: List of numbers that will be used as initial values for iteration.
* LINE4: A list of space-separated labels of each one of the nuclei
* LINE5 to the end: The pH (or pD) value of this point followed by the NMR shift
    of each one of the nuclei. If some data is missing, those numbers can be
    replaced with an **X**.

.. literalinclude:: input_example
   :caption: Input File Example
   :linenos:

.. seealso:: :mod:`gems.libio.load_file`

Interpreting the output
=======================

Bear in mind that :program:`GEMS` is a least-squares minimisation algorithm. In order to
have meaningful results you have to be careful with the data you feed the program
with and analyse carefuly the output provided. 

:program:`GEMS`  provided a verbose output with lots of information that will help you
to interpret what happened inside.

Refinement result
******************

After the optimization ends, the refined parameters are shown with the calculated errors.
They are displayed in order.

The correlation matrix shows how the refined parameters are related to each other.
A good model should show a low correlation between the parameters. If some 
parameters show high correlation it usually means that the model is oversized.

Microstates and Microsteps Analysis
***********************************

For each one of the macrostates, the possible microstates are shown. For each one
of the microstates, the refined free energy *F* with error is shown along with
the calculated population of that microstate.

Calculated Macroconstants
*************************

The macroconstants are displayed. These numbers represent the protonation
(or deuteration) equilibria. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
   tutorial
   algorithm
   api
   changelog

