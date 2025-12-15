"""
Set up simulations of two 1st order processes:
    - A -> B
    - A ->
    - A -> C + D
    - A -> 2C
Each process is simulated with three different population structures
of species A with respect to the rate constant `k`:
    - Homogeneous population
    - Heterogeneous population: two distinct subspecies
    - Heterogeneous population: normally-distributed k values
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

from abstochkin.base import AbStochKin
from abstochkin.utils import r_squared
from het_predictions_o1 import *


class TestFirstOrderSimulation1(unittest.TestCase):
    """ A -> B or A -> , homogeneous population """

    def setUp(self):
        """
        Homogeneous population structure of species A with respect to the
        rate constant `k`.
        """
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A -> B", 0.3)
        self.sim1a.simulate(p0={'A': 100, 'B': 0}, t_max=15, dt=0.1, show_plots=False,
                            max_agents_by_species={'A': 100, 'B': 100})

        self.sim1b = AbStochKin()
        self.sim1b.add_process_from_str("A -> ", 0.5)
        self.sim1b.simulate(p0={'A': 10}, t_max=15, dt=0.1, n=200, show_plots=False,
                            max_agents_by_species={'A': 10})
        
        self.sim1c = AbStochKin()
        self.sim1c.add_process_from_str("A -> B + C", 0.15)
        self.sim1c.simulate(p0={'A': 30, 'B': 0, 'C': 0}, t_max=15, dt=0.01, n=100,
                            show_plots=False, max_agents_by_species={'A': 30, 'B': 30, 'C': 30})
        
        self.sim1d = AbStochKin()
        self.sim1d.add_process_from_str("A -> 2C", 0.25)
        self.sim1d.simulate(p0={'A': 50, 'C': 0}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 50, 'C': 100})

    def test_simulations_1(self):
        """ Test simulations of process with homogeneous populations. """
        self.assertEqual(self.sim1a.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(self.sim1a.sims[0].k_vals['A -> B, k = 0.3']), {0.3})
        self.assertGreaterEqual(self.sim1a.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1a.sims[0].results['B']['R^2'], 0.999)

        self.assertEqual(self.sim1b.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(self.sim1b.sims[0].k_vals['A -> , k = 0.5']), {0.5})
        self.assertGreaterEqual(self.sim1b.sims[0].results['A']['R^2'], 0.9985)
        
        self.assertEqual(self.sim1c.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(self.sim1c.sims[0].k_vals['A -> B + C, k = 0.15']), {0.15})
        self.assertGreaterEqual(self.sim1c.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1c.sims[0].results['B']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1c.sims[0].results['C']['R^2'], 0.999)
        
        self.assertEqual(self.sim1d.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(self.sim1d.sims[0].k_vals['A -> 2 C, k = 0.25']), {0.25})
        self.assertGreaterEqual(self.sim1d.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1d.sims[0].results['C']['R^2'], 0.999)


class TestFirstOrderSimulation2(unittest.TestCase):
    """ A -> B or A -> , heterogeneous population - distinct subspecies """

    def setUp(self):
        """
        Heterogeneous population structure of species A with respect to the
        rate constant `k` with two distinct subspecies.
        """
        self.sim2a = AbStochKin()
        self.sim2a.add_process_from_str("A -> B", [0.1, 0.3])
        self.sim2a.simulate(p0={'A': 10, 'B': 0}, t_max=10, dt=0.1, n=200, show_plots=False,
                            max_agents_by_species={'A': 10, 'B': 10})

        self.sim2b = AbStochKin()
        self.sim2b.add_process_from_str("A -> ", [0.5, 0.2])
        self.sim2b.simulate(p0={'A': 16}, t_max=20, dt=0.1, n=200, show_plots=False,
                            max_agents_by_species={'A': 16})
        
    def test_simulations_2(self):
        """ Test simulations of process with heterogeneous population of A. """
        self.assertEqual(self.sim2a.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim2a.sims[0].k_het_metrics['A -> B, k = [0.1, 0.3]']), 7)
        self.assertSetEqual(set(self.sim2a.sims[0].k_vals['A -> B, k = [0.1, 0.3]']), {0.3, 0.1})
        # Now compare simulated average trajectory to theoretical expectation
        r_sq_avg2a = r_squared(self.sim2a.sims[0].results['A']['N_avg'], sim2a_theory_avg)
        self.assertGreaterEqual(r_sq_avg2a, 0.9975)
        r_sq_std2a = r_squared(self.sim2a.sims[0].results['A']['N_std'], sim2a_theory_std)
        self.assertGreaterEqual(r_sq_std2a, 0.88)

        self.assertEqual(self.sim2b.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim2b.sims[0].k_het_metrics['A -> , k = [0.5, 0.2]']), 7)
        self.assertSetEqual(set(self.sim2b.sims[0].k_vals['A -> , k = [0.5, 0.2]']), {0.2, 0.5})
        # Now compare simulated average trajectory to theoretical expectation
        r_sq_avg2b = r_squared(self.sim2b.sims[0].results['A']['N_avg'], sim2b_theory_avg)
        self.assertGreaterEqual(r_sq_avg2b, 0.999)
        r_sq_std2b = r_squared(self.sim2b.sims[0].results['A']['N_std'], sim2b_theory_std)
        self.assertGreaterEqual(r_sq_std2b, 0.9)


class TestFirstOrderSimulation3(unittest.TestCase):
    """ A -> B or A -> , heterogeneous population - normally-distributed `k`. """

    def setUp(self):
        """
        Heterogeneous population structure of species A with respect to the
        rate constant `k` with normally-distributed `k` values.
        """
        self.sim3a = AbStochKin()
        self.sim3a.add_process_from_str("A -> B", (0.3, 0.1))
        self.sim3a.simulate(p0={'A': 10, 'B': 0}, t_max=10, dt=0.1, n=200, show_plots=False,
                            max_agents_by_species={'A': 10, 'B': 10})

        self.sim3b = AbStochKin()
        self.sim3b.add_process_from_str("A -> ", (0.2, 0.05))
        self.sim3b.simulate(p0={'A': 15}, t_max=20, dt=0.1, show_plots=False,
                            max_agents_by_species={'A': 15})

    def test_simulations_3(self):
        """ Test simulations of process with normally-distributed k values. """
        self.assertEqual(self.sim3a.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim3a.sims[0].k_het_metrics['A -> B, k = (0.3, 0.1)']), 7)
        r_sq_avg3a = r_squared(self.sim3a.sims[0].results['A']['N_avg'], sim3a_theory_avg)
        self.assertGreaterEqual(r_sq_avg3a, 0.995)
        r_sq_std3a = r_squared(self.sim3a.sims[0].results['A']['N_std'], sim3a_theory_std)
        self.assertGreaterEqual(r_sq_std3a, 0.9)

        self.assertEqual(self.sim3b.sims[0]._het_processes_num, 1)
        self.assertEqual(len(self.sim3b.sims[0].k_het_metrics['A -> , k = (0.2, 0.05)']), 7)
        r_sq_avg3b = r_squared(self.sim3b.sims[0].results['A']['N_avg'], sim3b_theory_avg)
        self.assertGreaterEqual(r_sq_avg3b, 0.998)
        r_sq_std3b = r_squared(self.sim3b.sims[0].results['A']['N_std'], sim3b_theory_std)
        self.assertGreaterEqual(r_sq_std3b, 0.9)


if __name__ == '__main__':
    unittest.main()
