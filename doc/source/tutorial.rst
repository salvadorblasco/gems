Tutorials
#########

Tutorial 1. Using GEMS from command line
****************************************

.. note:: Updated for GEMS 0.4

If you are familiar with Linux, you must be familiar with the use of
a terminal prompt.

Command line options
====================

From the prompt you can type :command:`$ ./gems.py -h` will print all available
options.

.. literalinclude:: console_options

The most important option is :command:`-f FILE` which loads an input file
:command:`FILE`, starts the calculation and output the result.

.. _info_symmetry:

Tutorial 2: Obtaining information based on molecule symmetry
************************************************************

:program:`GEMS` when used with the option :samp:`-s` displays information about the microstates
and microsteps based on the symmetry provided. The symmetry of the molecule is described
by using characters to tag each one of the protonation centres and a number indicating
the multiplicity of said centres. For example, :samp:`A2B` indicates a molecule with three
protonation centres where two of them are chemically equivalent. Any letter can be used
to indicate a protonation centre. :samp:`A2B`, :samp:`X2P`, :samp:`BA2` and :samp:`AAB` 
are equivalent and should yield the same result except with different labels and in 
different order. ::

    molecule symmetry:  A2B 

    n = 0, 1 microstate(s):  Ø
       Ø  +  A  ->  A
       Ø  +  B  ->  B
      2 microstep(s)

    n = 1, 2 microstate(s):  B, A
       B  +  A  ->  AB
       A  +  A  ->  A2
       A  +  B  ->  AB
      3 microstep(s)

    n = 2, 2 microstate(s):  AB, A2
       AB  +  A  ->  A2B
       A2  +  B  ->  A2B
      2 microstep(s)

    n = 3, 1 microstate(s):  A2B

    total microstates = 6
    total microsteps = 7

    Fitting information
      first order parameters:   2
      second order parameters:  2
      third order parameters:   1
      - - - - - - - - - - - - - -
      total parameters:         5

The output is a list of each one of the macrostates, from *n* =0 to *n* =3. Inside each
one of the macrostates there is a list of the possible microstates and microsteps.

You can try with any symmetry.

.. code-block:: shell

    $./rungems.py -s A4B2
    $./rungems.py -s A6B2C
