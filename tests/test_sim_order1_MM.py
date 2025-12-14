"""
Set up simulations of a 1st order Michaelis-Menten process:
    - A -> B, catalyzed by species E
    - A -> , catalyzed by species C

Each process is simulated with three different population structures
of species A with respect to the rate constant `Km`:
    - Homogeneous population
    - Heterogeneous population: two distinct subspecies of `Km` values
    - Heterogeneous population: normally-distributed `Km` values
"""

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

from numpy import nan, mean, std, int64

from abstochkin.base import AbStochKin


class TestFirstOrderMMSimulation1(unittest.TestCase):
    def setUp(self):
        """
        Homogeneous population structure of species A with respect to the
        rate constant `k` and `Km`.
        """
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A -> B", 0.5, catalyst='E', Km=10)
        self.sim1a.simulate(p0={'A': 50, 'B': 0, 'E': 10}, t_max=20, dt=0.01, n=100,
                            show_plots=False, max_agents_by_species={'A': 50, 'B': 50, 'E': 10})

        self.sim1b = AbStochKin()
        self.sim1b.add_process_from_str("A -> ", 0.3, catalyst='C', Km=10)
        self.sim1b.simulate(p0={'A': 40, 'C': 10}, t_max=20, dt=0.01, n=100, show_plots=False)

    def test_simulations_1(self):
        """ Test simulations of process with homogeneous populations. """
        self.assertEqual(self.sim1a.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1a.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1a.sims[0].results['B']['R^2'], 0.999)
        self.assertIs(self.sim1a.sims[0].results['E']['R^2'], nan)

        self.assertEqual(self.sim1b.sims[0]._het_processes_num, 0)
        ''' Make sure max agents of catalyst (unmodified species) are not 
        determined by using the multiplier when `max_agents_by_species` is 
        not specified. (see below assertion) '''
        self.assertEqual(self.sim1b.sims[0].max_agents['C'], self.sim1b.sims[0].p0['C'])
        self.assertGreaterEqual(self.sim1b.sims[0].results['A']['R^2'], 0.999)
        self.assertIs(self.sim1b.sims[0].results['C']['R^2'], nan)


class TestFirstOrderMMSimulation2(unittest.TestCase):
    def setUp(self):
        """
        Heterogeneous population structure of species A with respect `Km`:
        two distinct subspecies or normal distribution.
        """
        self.sim2a = AbStochKin()
        self.sim2a.add_process_from_str("A -> B", (0.25, 0.05),
                                        catalyst='E', Km=[10, 15])
        self.sim2a.simulate(p0={'A': 50, 'B': 0, 'E': 10}, t_max=20, run=False)
        self.sim2a.sims[0]._setup_data()

        self.sim2b = AbStochKin()
        self.sim2b.add_process_from_str("A -> B", [0.25, 0.15],
                                        catalyst='E', Km=(10, 2))
        self.sim2b.simulate(p0={'A': 50, 'B': 0, 'E': 10}, t_max=20, run=False)
        self.sim2b.sims[0]._setup_data()

    def test_simulations_2(self):
        """ Test simulations of process with heterogeneous populations. """
        proc2a = self.sim2a.processes[0]
        self.assertEqual(self.sim2a.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim2a.sims[0].Km_vals), 1)
        self.assertEqual(mean(self.sim2a.sims[0].Km_vals[proc2a]),
                         mean(self.sim2a.processes[0].Km))
        self.assertEqual(len(self.sim2a.sims[0].Km_het_metrics), 1)

        proc2b = self.sim2b.processes[0]
        self.assertEqual(self.sim2b.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim2b.sims[0].Km_vals), 1)
        self.assertEqual(self.sim2b.sims[0].Km_vals[proc2b].dtype, int64)
        avg_vals = mean(self.sim2b.sims[0].Km_vals[proc2b])
        self.assertLessEqual(abs(avg_vals - self.sim2b.processes[0].Km[0]), 1)
        std_vals = std(self.sim2b.sims[0].Km_vals[proc2b])
        self.assertLessEqual(abs(std_vals - self.sim2b.processes[0].Km[1]), 0.25)


if __name__ == '__main__':
    unittest.main()
