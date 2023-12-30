#  Copyright (c) 2024, Alex Plakantonakis.
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

import unittest
import numpy as np
from abstochkin.utils import rng_streams, macro_to_micro


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


if __name__ == '__main__':
    unittest.main()
