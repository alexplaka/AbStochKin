#  Copyright (c) 2024-2026, Alex Plakantonakis.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import io
import logging
import unittest
from contextlib import redirect_stdout
from time import sleep

import numpy as np

from abstochkin.utils import rng_streams, macro_to_micro, measure_runtime, r_squared, log_exceptions


class TestRng(unittest.TestCase):
    """
    Make sure random number streams spawned from the same seed
    are identical. Also, make sure random number streams spawned
    from different seeds are different.
    """

    def setUp(self):
        self.n_gens = 10  # number of generators/streams
        self.streams_1 = rng_streams(self.n_gens, random_state=19)
        self.streams_2 = rng_streams(self.n_gens, random_state=19)
        self.streams_3 = rng_streams(self.n_gens, random_state=42)
        self.num = 100  # number of random numbers for each generator/stream to generate

    def test_rng_streams_integers(self):
        test_1 = np.zeros([self.n_gens, self.num])
        test_2 = np.zeros([self.n_gens, self.num])
        test_3 = np.zeros([self.n_gens, self.num])
        for i in range(self.n_gens):
            test_1[i, :] = self.streams_1[i].integers(0, 1000, self.num)
            test_2[i, :] = self.streams_2[i].integers(0, 1000, self.num)
            test_3[i, :] = self.streams_3[i].integers(0, 1000, self.num)
        self.assertEqual(np.sum(test_1 - test_2), 0)
        self.assertNotEqual(np.sum(test_1 - test_3), 0)

    def test_rng_streams_floats(self):
        test_1 = np.zeros([self.n_gens, self.num])
        test_2 = np.zeros([self.n_gens, self.num])
        test_3 = np.zeros([self.n_gens, self.num])
        for i in range(self.n_gens):
            test_1[i, :] = self.streams_1[i].random(self.num)
            test_2[i, :] = self.streams_2[i].random(self.num)
            test_3[i, :] = self.streams_3[i].random(self.num)
        self.assertEqual(np.sum(test_1 - test_2), 0)
        self.assertNotEqual(np.sum(test_1 - test_3), 0)

    def test_rng_streams_normal(self):
        test_1 = self.streams_1[0].normal(0, 1, self.num)
        test_2 = self.streams_2[0].normal(0, 1, self.num)
        test_3 = self.streams_3[0].normal(0, 1, self.num)
        self.assertEqual(np.sum(test_1 - test_2), 0)
        self.assertNotEqual(np.sum(test_1 - test_3), 0)

    def test_rng_streams_edge_cases(self):
        # Test n=0
        streams_0 = rng_streams(0, random_state=19)
        self.assertEqual(len(streams_0), 0)

        # Test n=1
        streams_1 = rng_streams(1, random_state=19)
        self.assertEqual(len(streams_1), 1)
        self.assertIsInstance(streams_1[0], np.random.Generator)


class TestMacroToMicro(unittest.TestCase):
    def test_conversion(self):
        self.assertAlmostEqual(macro_to_micro(0.001, 1e-6, 0),
                               6.02214076e14)
        self.assertAlmostEqual(macro_to_micro(1e-9, 1e-6, 0),
                               6.02214076e8)

        self.assertEqual(macro_to_micro(0.1, 1e-6, 1),
                         0.1)
        self.assertEqual(macro_to_micro(0.05, 1e-6, 1),
                         0.05)

        self.assertAlmostEqual(macro_to_micro(10, 1e-8, 2),
                               1.660539e-15)
        self.assertAlmostEqual(macro_to_micro(0.01, 1e-15, 2),
                               1.660539e-12)

    def test_vectorized_conversion(self):
        result_0 = macro_to_micro([0.001, 0.002, 0.003], 1e-6, 0)
        expected_0 = [6.02214076e14, 1.204428152e+15, 1.806642228e+15]
        for i in range(len(expected_0)):
            self.assertAlmostEqual(result_0[i], expected_0[i], places=6)

        result_1 = macro_to_micro([0.15, 0.11], 2e-5, 1)
        expected_1 = [0.15, 0.11]
        self.assertListEqual(result_1, expected_1)

        result_2 = macro_to_micro((0.001, 0.00025), 1e-10, 2)
        expected_2 = (1.6605390671738466e-17, 4.1513476679346165e-18)
        for i in range(len(expected_2)):
            self.assertAlmostEqual(result_2[i], expected_2[i], places=6)

    def test_inverse_conversion(self):
        self.assertEqual(macro_to_micro(60221.4076, 1e-16, 0, inverse=True), 0.001)
        self.assertEqual(macro_to_micro(0.1, 1e-6, 1, inverse=True), 0.1)
        self.assertEqual(macro_to_micro(8.302695335869234e-05, 1e-20, 2, inverse=True), 0.5)

        # Test with analogy to composition of functions $f(f^{-1}(x))=x$
        k_micro_0 = 1.5  # 0th order microscopic k value
        self.assertEqual(
            macro_to_micro(
                macro_to_micro(k_micro_0, volume=1e-10, order=0, inverse=True),  # get macroscopic value
                volume=1e-10, order=0),
            k_micro_0)

        k_micro_1 = 0.5  # 1st order microscopic k value
        self.assertEqual(
            macro_to_micro(
                macro_to_micro(k_micro_1, volume=1e-5, order=1, inverse=True),  # get macroscopic value
                volume=1e-5, order=1),
            k_micro_1)

        k_micro_2 = 0.002  # 2nd order microscopic k value
        self.assertEqual(
            macro_to_micro(
                macro_to_micro(k_micro_2, volume=1e-15, order=2, inverse=True),  # get macroscopic value
                volume=1e-15, order=2),
            k_micro_2)


class TestMeasureRuntime(unittest.TestCase):
    def test_measure_runtime_seconds(self):
        @measure_runtime
        def sleeping_function():
            sleep(1.243)

        output = io.StringIO()
        handler = logging.StreamHandler(output)
        logger = logging.getLogger()
        logger.addHandler(handler)

        with redirect_stdout(output):
            sleeping_function()

        self.assertEqual(output.getvalue(), "Simulation Runtime: 1.243 sec\n")


class TestRSquared(unittest.TestCase):
    def test_perfect_fit(self):
        actual = np.array([1.0, 2.0, 3.0, 4.0])
        theoretical = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertEqual(r_squared(actual, theoretical), 1.0)

    def test_no_fit(self):
        actual = np.array([1.0, 1.0, 1.0, 1.0])
        theoretical = np.array([2.0, 2.0, 2.0, 2.0])
        self.assertTrue(np.isnan(r_squared(actual, theoretical)))

    def test_partial_fit(self):
        actual = np.array([1.0, 2.0, 3.0, 4.0])
        theoretical = np.array([1.5, 2.5, 3.5, 4.5])
        r2 = r_squared(actual, theoretical)
        self.assertTrue(0 < r2 < 1)

    def test_with_nan_values(self):
        actual = np.array([1.0, np.nan, 3.0, 4.0])
        theoretical = np.array([1.0, 2.0, 3.0, 4.0])
        r2 = r_squared(actual, theoretical)
        self.assertFalse(np.isnan(r2))

    def test_different_lengths(self):
        actual = np.array([1.0, 2.0, 3.0])
        theoretical = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            r_squared(actual, theoretical)

    def test_all_nan_values(self):
        actual = np.array([np.nan, np.nan, np.nan])
        theoretical = np.array([1.0, 2.0, 3.0])
        r2 = r_squared(actual, theoretical)
        self.assertTrue(np.isnan(r2))

    def test_empty_arrays(self):
        actual = np.array([])
        theoretical = np.array([])
        r2 = r_squared(actual, theoretical)
        self.assertTrue(np.isnan(r2))


class TestMacroToMicroErrors(unittest.TestCase):
    def test_negative_volume(self):
        with self.assertRaises(AssertionError):
            macro_to_micro(0.001, -1e-6, 0)

    def test_negative_order(self):
        with self.assertRaises(AssertionError):
            macro_to_micro(0.001, 1e-6, -1)

    def test_negative_values_list(self):
        with self.assertRaises(AssertionError):
            macro_to_micro([-0.001, 0.002], 1e-6, 0)

    def test_mixed_types(self):
        result = macro_to_micro([1, 2.5], 1e-6, 0)
        self.assertIsInstance(result[0], float)
        self.assertIsInstance(result[1], float)

    def test_empty_sequence(self):
        result_list = macro_to_micro([], 1e-6, 0)
        self.assertEqual(result_list, [])
        result_tuple = macro_to_micro((), 1e-6, 0)
        self.assertEqual(result_tuple, ())

    def test_zero_values(self):
        self.assertEqual(macro_to_micro(0, 1e-6, 0), 0)
        self.assertEqual(macro_to_micro([0, 0], 1e-6, 1), [0, 0])
        self.assertEqual(macro_to_micro((0, 0), 1e-6, 2), (0, 0))


class TestLogExceptions(unittest.TestCase):
    def test_no_exception(self):
        output = io.StringIO()
        handler = logging.StreamHandler(output)
        logger = logging.getLogger()
        logger.addHandler(handler)

        with log_exceptions():
            x = 1 + 1

        self.assertEqual(output.getvalue(), "")

    def test_with_exception(self):
        output = io.StringIO()
        handler = logging.StreamHandler(output)
        logger = logging.getLogger()
        logger.addHandler(handler)

        with self.assertRaises(ZeroDivisionError):
            with log_exceptions():
                x = 1 / 0

        self.assertIn("ZeroDivisionError", output.getvalue())


if __name__ == '__main__':
    unittest.main()
