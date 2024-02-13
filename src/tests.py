#!/usr/bin/env python3

# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
# ! This file is part of the GEMS software                                   !
# !                                                                          !
# ! Copyright (c) 2020-2024 by Salvador Blasco <salvador.blasco@gmail.com>   !
# ! Licensed under MIT license (see file LICENSE)                            !
# +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

"""+-------------------------------+
   |    Testing suite for GEMS     |
   +-------------------------------+ """

import unittest
import unittest.mock

import contextlib
import io
import sys

import numpy as np


@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.StringIO()
    yield
    sys.stdout = save_stdout


class TestLibFit(unittest.TestCase):
    'Unit test for functions in fit.py'

    @unittest.skip('broken')
    def test_fit_shifts(self):
        from gems import fit

        with np.load('testsdata.npz') as f:
            shifts = f['yvalues']
            macrop = f['macrostate_probability']

        calc_shifts = fit.fit_shifts(macrop, shifts)        # <= TEST CALL
        np.testing.assert_array_almost_equal(calc_shifts, shifts, decimal=0)

        shifts2 = np.ma.array(shifts)
        for r, c in ((1, 2), (5, 6), (0, 5)):
            shifts2[r, c] = np.ma.masked
        calc_shifts = fit.fit_shifts(macrop, shifts)        # <= TEST CALL
        np.testing.assert_array_almost_equal(calc_shifts, shifts, decimal=0)

    @unittest.skip('broken')
    def test_fit_shifts2(self):
        from gems import fit
        ftest = fit.fit_shifts2

        with np.load('testsdata.npz') as f:
            shifts = f['yvalues']
            theta = f['theta']

        calc_shifts = ftest(theta, shifts)        # <= TEST CALL
        np.testing.assert_array_almost_equal(calc_shifts, shifts, decimal=2)

        shifts2 = np.ma.array(shifts)
        for r, c in ((1, 2), (5, 6), (0, 5)):
            shifts2[r, c] = np.ma.masked
        calc_shifts = ftest(theta, shifts)       # <= TEST CALL
        np.testing.assert_array_almost_equal(calc_shifts, shifts, decimal=0)

    @unittest.skip('broken')
    def test_postfit(self):
        import gems.libio
        result = gems.libio.import_numpy('testsdata.npz')
        indict = {k: result[k] for k in ('result', 'molecule', 'xvalues', 'yvalues', 'labels', 'order')}
        
        import gems.fit
        gems.fit.postfit(indict)
        for k in result:
            if k == 'microstate_probability':       # TODO not implemented, fails comparing numpy arrays
                continue
            with self.subTest(k=k):
                if isinstance(result[k], np.ndarray):
                    np.testing.assert_array_almost_equal(result[k], indict[k])
                elif isinstance(result[k], dict):
                    self.assertDictEqual(result[k], indict[k])
                else:
                    self.assertEqual(result[k], indict[k])


class TestLibUk(unittest.TestCase):
    def setUp(self):
        self.energy_sample = {(0, 0, 0): 0.0, (1, 0, 0): -20.0,
                              (0, 1, 0): -20.0, (0, 0, 1): -25.0,
                              (1, 1, 0): -40.0, (1, 0, 1): -45.0,
                              (0, 1, 1): -45.0, (1, 1, 1): -60.0}

    def test_2_equal(self):

        def two_equal(seq):
            import collections
            count = collections.Counter(seq)
            return any(i > 1 for i in count.values())

        seq_true = ((0, 0), (1, 1), (0, 1, 1), (2, 2, 1))
        for seq in seq_true:
            with self.subTest(seq=seq):
                self.assertTrue(two_equal(seq))

        seq_false = ((2, 0), (1, 2), (1, 0, 3), (3, 2, 1))
        for seq in seq_false:
            with self.subTest(seq=seq):
                self.assertFalse(two_equal(seq))


    def test_avg_prtdeg(self):
        import gems.libuk
        
        mstp = {(0, 0, 0): 8.756505387377377e-27,
                (1, 0, 0): 4.248351647373827e-18,
                (0, 1, 0): 4.248351647373827e-18,
                (0, 0, 1): 6.305112889653852e-16,
                (1, 1, 0): 2.061152357167627e-09,
                (1, 0, 1): 3.0590213271896984e-07,
                (0, 1, 1): 3.0590213271896984e-07,
                (1, 1, 1): 0.9999993861345816}

        reslt1 = gems.libuk.avg_prtdeg(mstp)
        reslt2 = np.array([[0.33333323], [0.33333323], [0.33333333]])
        np.testing.assert_array_almost_equal(reslt1, reslt2)

    @unittest.skip("broken")
    def test_compute_mapping(self):
        import gems.libuk
        test_molec = "A2B2"
        mapping = gems.libuk.compute_mapping(test_molec, isomorphisms)

    @unittest.skip("broken")
    def test_microstate_multiplicity(self):
        data = (('AABC', 'AB', 2),
                ('AABC', 'AA', 1),
                ('AABB', 'AB', 4),
                ('AABB', 'BA', 4),
                ('AABC', 'AC', 2))
        from gems import libuk
        for molec, micrst, result in data:
            with self.subTest(molec=molec, microstate=micrst):
                out = libuk.microstate_multiplicity(molec, micrst)
                self.assertEqual(result, out)

    @unittest.skip("broken")
    def test_order_terms(self):
        import gems.libuk
        self.assertListEqual(gems.libuk.order_terms('AABC', 2),
                             ['AA', 'AB', 'AC', 'BC'])

    def test_conditional_probability2(self):
        from math import exp
        from gems import libuk
        energy = self.energy_sample
        consts = [1.0, 2*exp(20.0) + exp(25.0), exp(40.0) + 2*exp(45.0),
                  exp(60.0)]
        calc_condprob = {k: exp(-v)/consts[sum(k)] for k, v in energy.items()}
        test_condprob = libuk.conditional_probability2(energy)
        self.assertDictEqual(calc_condprob, test_condprob)

    def test_kdeflation(self):
        from gems import libuk
        n = range(-1, 10)
        y = ((), (), (1.0,), (1.0, 0.25),
             (1.0, 0.3333333333333333, 0.1111111111111111),
             (1.0, 0.375, 0.16666666666666666, 0.0625),
             (1.0, 0.4, 0.2, 0.1, 0.04),
             (1.0, 0.4166666666666667, 0.2222222222222222, 0.125,
              0.06666666666666667, 0.027777777777777776),
             (1.0, 0.42857142857142855, 0.2380952380952381,
              0.14285714285714285, 0.08571428571428572, 0.047619047619047616,
              0.02040816326530612),
             (1.0, 0.4375, 0.25, 0.15625, 0.1, 0.0625, 0.03571428571428571,
              0.015625),
             (1.0, 0.4444444444444444, 0.2592592592592593, 0.16666666666666666,
              0.1111111111111111, 0.07407407407407407, 0.047619047619047616,
              0.027777777777777776, 0.012345679012345678))
        for n, r in zip(n, y):
            for a, b in zip(r, libuk.kdeflation(n)):
                self.assertEqual(a, b)

    def test_macro_constants(self):
        from math import exp, log
        from gems import libuk
        energy = self.energy_sample
        test_consts = [1.0, 2*exp(20.0) + exp(25.0), exp(40.0) + 2*exp(45.0),
                       exp(60.0)]
        calc_consts = libuk.macro_constants(energy)
        for n, (a, b) in enumerate(zip(test_consts, calc_consts)):
            with self.subTest(n=n):
                self.assertAlmostEqual(log(a), log(b), places=7)

    def test_matrix_a(self):
        from gems import libuk
        condprob = {(0, 0, 0): 1.0,
                    (1, 0, 0): 0.006648354478866004,
                    (0, 1, 0): 0.006648354478866004,
                    (0, 0, 1): 0.986703291042268,
                    (1, 1, 0): 0.003357661626502615,
                    (1, 0, 1): 0.4983211691867487,
                    (0, 1, 1): 0.4983211691867487,
                    (1, 1, 1): 1.0}
        calc_a = np.array(((0.0, 0.006648354478866004, 0.003357661626502615 +
                            0.4983211691867487, 1.0),
                           (0.0, 0.006648354478866004, 0.003357661626502615 +
                            0.4983211691867487, 1.0),
                           (0.0, 0.986703291042268, 2*0.4983211691867487,
                            1.0)))
        test_a = libuk.matrix_a(condprob, 3)
        np.testing.assert_array_almost_equal(calc_a, test_a)

    @unittest.skip("broken")
    def test_name_microstate(self):
        from gems import libuk
        feed = (('AABC', (1, 0, 0, 1), 'AC'),
                ('AAB', (1, 0, 0), 'A'),
                ('AABB', (1, 0, 1, 1), 'ABB'))
        for molec, ms, result in feed:
            retv = libuk.name_microstate(molec, ms)
            self.assertEqual(result, retv)

    def test_name_terms(self):
        from gems import libuk
        feed = (('AABC', 1, ['A', 'B', 'C']),
                ('AABC', 2, ['AA', 'AB', 'AC', 'BC']),
                ('AABC', 3, ['AAB', 'AAC', 'ABC']))
        for molec, lvl, result in feed:
            retv = libuk.name_terms(molec, lvl)
            self.assertEqual(result, retv)

    def test_error_macro_constants(self):
        from math import exp
        from gems import libuk
        energy = self.energy_sample
        errenergy = {(0, 0, 0): 0.0, (1, 0, 0): 0.1, (0, 1, 0): 0.1,
                     (0, 0, 1): 0.1, (1, 1, 0): 0.1, (1, 0, 1): 0.1,
                     (0, 1, 1): 0.1, (1, 1, 1): 0.1}
        # test_consts = [1.0, 2*exp(20.0) + exp(25.0), exp(40.0) + 2*exp(45.0),
        #                exp(60.0)]
        test_ecnsts = [0.0, .2*exp(2*20.0) + .1*exp(2*25.0), .1*exp(2*40.0) +
                       .2*exp(2*45.0), .1*exp(2*60.0)]
        calc_econsts = libuk.error_macro_constants(energy, errenergy)
        for n, (a, b) in enumerate(zip(test_ecnsts[1:], calc_econsts[1:])):
            with self.subTest(n=n):
                self.assertAlmostEqual(a/b, 1.0, places=7)

    def test_microstate_probability(self):
        from math import exp
        from gems import libuk
        energy = self.energy_sample
        ah = 1.0
        msts = {mv: exp(-e) for mv, e in energy.items()}
        norm = sum(msts.values())
        calc_microp = libuk.microstate_probability(energy, ah)
        test_microp = {ms: e/norm for ms, e in msts.items()}
        for n, (ms, p) in enumerate(test_microp.items()):
            with self.subTest(micros_state=ms):
                self.assertAlmostEqual(p, calc_microp[ms], places=7)

    def test_expand_name(self):
        from gems import libuk
        data_in = ('A2', 'AB', 'A3', 'A2B', 'ABC', 'A4', 'A3B', 'A2B2', 'A2BC',
                   'ABCD', 'A5', 'A4B', 'ABCDE', 'A6', 'A5B', 'A4B2', 'ABCDEF')
        data_out = ('AA', 'AB', 'AAA', 'AAB', 'ABC', 'AAAA', 'AAAB', 'AABB',
                    'AABC', 'ABCD', 'AAAAA', 'AAAAB', 'ABCDE', 'AAAAAA',
                    'AAAAAB', 'AAAABB', 'ABCDEF')
        for i, o in zip(data_in, data_out):
            assert o == libuk.expand_name(i)

    def test_num_prot_centres(self):
        from gems import libuk
        energy = self.energy_sample
        n = libuk.num_prot_centres(energy)
        self.assertEqual(n, 3)


class TestIO(unittest.TestCase):
    """Set of tests for :file:`libio.py`."""

    def test_export_numpy(self):
        infodict = {'a': 1, 'b': 2, 'c': 3}
        import gems.libio
        with unittest.mock.patch('numpy.savez_compressed') as mock_savez:
            gems.libio.export_numpy('example.npz', infodict)
            mock_savez.assert_called_with('example.npz', **infodict)

    def test_test_modules(self):
        import gems.libio
        with nostdout():
            gems.libio.test_modules()

    def test_print_welcome(self):
        import gems.libio
        with nostdout():
            gems.libio.print_welcome()

    @unittest.skip('not implemented')
    def test_load_file_spectra(self):
        pass

    @unittest.skip('not implemented')
    def test_load_file(self):
        pass


class TestIsomorph(unittest.TestCase):
    def test_clasify_microstates(self):
        import numpy as np
        import gems.isomorph
        isomorphisms = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                        np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])]
        microstates = {(1, 0, 0), (0, 0, 1), (0, 1, 0)}
        set1 = gems.isomorph.clasify_microstates(microstates, isomorphisms)
        set2 = {frozenset([(1, 0, 0)]), frozenset([(0, 0, 1), (0, 1, 0)])}
        self.assertSetEqual(set1, set2)

    def test_isomorphisms(self):
        import numpy as np
        import gems.isomorph
        ident = np.array([[0,1,1],[1,0,0],[1,0,0]], dtype=int)
        isom1 = gems.isomorph.find_isomorphisms(ident, 'BAA')
        isom2 = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                 np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])]
        for a, b in zip(isom1, isom2):
            np.testing.assert_array_equal(a, b)


class TestFittingParrams(unittest.TestCase):
    def test_class(self):
        import gems.vmc
        # parm = gems.vmc.FittingParams()

if __name__ == '__main__':
    unittest.main()
