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

import unittest
from abstochkin.het_calcs import get_het_processes, richness, idx_het, info_het
from abstochkin.process import Process, MichaelisMentenProcess, RegulatedProcess, ReversibleProcess
from numpy import array, log


class TestHeterogeneityCalculations(unittest.TestCase):

    def test_count_het_procs(self):
        a = [
            Process.from_string("X -> Y", k=(0.1, 0.02)),  # het
            Process.from_string("A -> B", k=0.5),  # hom
            Process.from_string("P -> ", k=[0.1, 0.15, 0.2]),  # het
            ReversibleProcess.from_string("C <-> D", k=0.1, k_rev=0.5),  # hom
            ReversibleProcess.from_string("C <-> D", k=[0.1, 0.15, 0.2], k_rev=0.5),  # het
            ReversibleProcess.from_string("C <-> D", k=(0.2, 0.1), k_rev=[0.5, 0.25]),  # het
            ReversibleProcess.from_string("C <-> D", k=[0.1, 0.15], k_rev=(0.5, 0.1)),  # het
            MichaelisMentenProcess.from_string("A -> B", k=0.1,
                                               catalyst='C', Km=10),  # hom
            MichaelisMentenProcess.from_string("A -> B", k=[0.1, 0.2],
                                               catalyst='C', Km=(20, 5)),  # het
            RegulatedProcess.from_string("A -> B", k=0.1,
                                         regulating_species='A', alpha=0,
                                         K50=10, nH=1),  # hom
            RegulatedProcess.from_string("A -> B", k=0.1,
                                         regulating_species='A', alpha=3,
                                         K50=(10, 2), nH=2),  # het
            RegulatedProcess.from_string("A -> C", k=0.3,
                                         regulating_species='A, C', alpha=[0, 2],
                                         K50=[15, 5], nH=[1, 2]),  # hom
            RegulatedProcess.from_string("A -> C", k=(0.3, 0.1),
                                         regulating_species='A, C', alpha=[0, 2],
                                         K50=[[15, 10], 5], nH=[1, 2])  # het
        ]
        self.assertEqual(len(get_het_processes(a)), 8)

    def test_richness(self):
        a = array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        b = array([0.1, 0.1, 0.10, 0.1, 0.100, 0.1, 0.1, 0.1, 0.25])
        c = array([1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3.0, 3, 3, 3])
        self.assertEqual(richness([]), 0)
        self.assertEqual(richness([0.1]), 1)
        self.assertEqual(richness(a), 1)
        self.assertEqual(richness(b), 2)
        self.assertEqual(richness(c), 3)

    def test_idx_het(self):
        self.assertEqual(idx_het([]), 0)
        self.assertEqual(idx_het([0]), 0)
        self.assertEqual(idx_het([1]), 0)
        self.assertEqual(idx_het([1, 1, 1, 1, 1, 1, 1]), 0)
        self.assertNotEqual(idx_het([1, 1, 1, 1, 1, 1, 1, 2]), 0)
        self.assertEqual(idx_het(array([1, 2])), 1)
        self.assertNotEqual(idx_het(array([1, 2, 3, 4, 5, 6, 6])), 1)
        self.assertAlmostEqual(idx_het(array([0.1] * 5 + [0.3] * 5)),
                               0.5555555555555556)
        self.assertAlmostEqual(idx_het(array([0.1] * 10000 + [0.3] * 10000)),
                               0.5000, places=4)
        self.assertEqual(idx_het(array([1, 2, 3, 4, 5, 6])), 1)

    def test_info_het(self):
        self.assertEqual(idx_het([]), 0)
        self.assertEqual(idx_het([0]), 0)
        self.assertEqual(idx_het([1]), 0)
        self.assertEqual(idx_het([1, 1, 1, 1, 1, 1, 1]), 0)
        self.assertNotEqual(idx_het([1, 1, 1, 1, 1, 1, 1, 2]), 0)
        self.assertEqual(info_het(array([1, 2])), log(2))
        self.assertEqual(info_het(array([1, 1, 2, 2])), log(2))
        self.assertEqual(info_het([1] * 10 + [2] * 10), log(2))
        self.assertGreater(info_het(array([1, 2, 3, 4, 5, 6, 7])),
                           info_het(array([1, 2, 3, 4, 5, 6])))


if __name__ == '__main__':
    unittest.main()
