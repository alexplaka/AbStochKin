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

import numpy as np

from abstochkin._simulation_methods import SimulationMethodsMixin
from abstochkin.base import AbStochKin


class TestSimulationMethodsMixin(unittest.TestCase):
    def setUp(self):
        self.func = SimulationMethodsMixin._o2_get_unique_pairs

        self.input_arrays = [
            np.array([(i, 0) for i in range(7)]),
            np.array([(i, 1) for i in range(7)]),
            np.concatenate((np.array([(i, 0) for i in range(7)]),
                            np.array([(i, 1) for i in range(7)]))),
            np.array([(0, 0), (0, 3), (1, 0), (1, 3), (1, 4)]),
            np.concatenate((np.array([(i, 0) for i in range(4)]),
                            np.array([(i, 1) for i in range(4)]),
                            np.array([(i, 2) for i in range(4)]),
                            np.array([(i, 3) for i in range(4)]))),
            np.array([(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]),
            np.array([(0, 0), (1, 1), (0, 2), (1, 3), (0, 4), (1, 5)])
        ]

        self.output = [self.func(arr) for arr in self.input_arrays]

        self.answer = [
            {(0, 0)},
            {(0, 1)},
            {(0, 0), (1, 1)},
            {(0, 0), (1, 3)},
            {(0, 0), (1, 1), (2, 2), (3, 3)},
            {(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)},
            {(0, 0), (1, 1)}
        ]

    def test__o2_get_unique_pairs(self):
        for out, ans in zip(self.output, self.answer):
            self.assertSetEqual(set([tuple(t) for t in out.tolist()]), ans)


class TestRegulatedProcessSetupData(unittest.TestCase):
    def setUp(self):
        # Only 1 regulating species
        self.sim0 = AbStochKin()
        self.sim0.add_process_from_str("A -> B", k=0.2,
                                       regulating_species='A', alpha=0,
                                       K50=15, nH=2)
        self.sim0.simulate(p0={'A': 50, 'B': 0}, t_max=10, run=False)
        self.sim0.sims[0]._setup_data()

        # 2 regulating species with homogeneity in K50
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A -> B", k=0.2,
                                        regulating_species=['A', 'B'], alpha=[0, 2],
                                        K50=[5, 10], nH=[2, 1])
        self.sim1a.simulate(p0={'A': 50, 'B': 0}, t_max=10, run=False)
        self.sim1a.sims[0]._setup_data()

        # 2 regulating species, both with heterogeneity in K50
        self.sim2a = AbStochKin()
        self.sim2a.add_process_from_str("A + B -> C", k=0.3,
                                        regulating_species='A, C', alpha=[2, 3],
                                        K50=[(10, 2), [15, 10]], nH=[1, 2])
        self.sim2a.simulate(p0={'A': 50, 'B': 50, 'C': 0}, t_max=10, run=False)
        self.sim2a.sims[0]._setup_data()

        # 2 regulating species, only one with heterogeneity in K50
        self.sim2b = AbStochKin()
        self.sim2b.add_process_from_str("A + B -> C", k=0.3,
                                        regulating_species='A, C', alpha=[2, 3],
                                        K50=[5, [15, 10]], nH=[1, 2])
        self.sim2b.simulate(p0={'A': 50, 'B': 50, 'C': 0}, t_max=10, run=False)
        self.sim2b.sims[0]._setup_data()

    def test_K50_vals_het_metrics(self):
        K50_vals_0 = self.sim0.sims[0].K50_vals[self.sim0.processes[0]]
        self.assertIsInstance(K50_vals_0, np.ndarray)
        self.assertEqual(len(K50_vals_0), 100)

        K50_vals_1a = self.sim1a.sims[0].K50_vals[self.sim1a.processes[0]]
        self.assertIsInstance(K50_vals_1a, list)
        self.assertEqual(len(K50_vals_1a), 2)
        self.assertEqual(np.mean(K50_vals_1a[0]), self.sim1a.processes[0].K50[0])
        self.assertEqual(np.std(K50_vals_1a[0]), 0)
        self.assertEqual(np.mean(K50_vals_1a[1]), self.sim1a.processes[0].K50[1])
        self.assertEqual(np.std(K50_vals_1a[1]), 0)
        K50_het_1a = self.sim1a.sims[0].K50_het_metrics[self.sim1a.processes[0]]
        self.assertIsInstance(K50_het_1a, list)
        self.assertEqual(len(K50_het_1a), 2)
        self.assertListEqual(K50_het_1a, [None, None])

        K50_vals_2a = self.sim2a.sims[0].K50_vals[self.sim2a.processes[0]]
        self.assertIsInstance(K50_vals_2a, list)
        self.assertEqual(len(K50_vals_2a), 2)
        self.assertLessEqual(abs(np.mean(K50_vals_2a[0]) - self.sim2a.processes[0].K50[0][0]), 1)
        self.assertLessEqual(abs(np.std(K50_vals_2a[0]) - self.sim2a.processes[0].K50[0][1]), 0.5)
        self.assertLessEqual(abs(np.mean(K50_vals_2a[1]) - np.mean(self.sim2a.processes[0].K50[1])), 1)
        self.assertLessEqual(np.std(K50_vals_2a[1]), 2.5)
        K50_het_2a = self.sim2a.sims[0].K50_het_metrics[self.sim2a.processes[0]]
        self.assertIsInstance(K50_het_2a, list)
        self.assertEqual(len(K50_het_2a), 2)
        self.assertNotEqual(K50_het_2a[0], None)
        self.assertNotEqual(K50_het_2a[1], None)

        K50_vals_2b = self.sim2b.sims[0].K50_vals[self.sim2b.processes[0]]
        self.assertIsInstance(K50_vals_2b, list)
        self.assertEqual(len(K50_vals_2b), 2)
        self.assertEqual(np.mean(K50_vals_2b[0]), self.sim2b.processes[0].K50[0])
        self.assertEqual(np.std(K50_vals_2b[0]), 0)
        self.assertLessEqual(abs(np.mean(K50_vals_2b[1]) - np.mean(self.sim2b.processes[0].K50[1])), 1)
        self.assertLessEqual(np.std(K50_vals_2b[1]), 2.5)
        K50_het_2b = self.sim2b.sims[0].K50_het_metrics[self.sim2b.processes[0]]
        self.assertIsInstance(K50_het_2b, list)
        self.assertEqual(len(K50_het_2b), 2)
        self.assertEqual(K50_het_2b[0], None)
        self.assertNotEqual(K50_het_2b[1], None)


if __name__ == '__main__':
    unittest.main()
