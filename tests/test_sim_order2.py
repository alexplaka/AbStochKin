"""
Set up simulations of 2nd order processes:
    - A + B -> C
    - A + A -> C
Each process is simulated with three different classes of 
inter-agent interactions wrt rate constant `k`:
    - Homogeneous population
    - Heterogeneous population: two k values (or interactions)
    - Heterogeneous population: normally-distributed k values
"""

import unittest

from numpy import unique, diagonal, sum

from abstochkin.base import AbStochKin


class TestSecondOrderSimulation1(unittest.TestCase):
    """ A + B -> C, homogeneous population
        A + B -> C + D, homogeneous population """

    def setUp(self):
        """
        Homogeneous population structure of interactions between species A
        and species B with respect to the rate constant `k`.
        """
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A + B -> C", k=0.03)
        self.sim1a.simulate(p0={'A': 50, 'B': 50, 'C': 0}, t_max=10, dt=0.01, n=100,
                            show_plots=False, max_agents_by_species={'A': 50, 'B': 50, 'C': 50})

        self.sim1b = AbStochKin()
        self.sim1b.add_process_from_str("A + B -> 2C", k=0.003)
        self.sim1b.simulate(p0={'A': 100, 'B': 70, 'C': 0}, t_max=10, dt=0.01, n=100,
                            show_plots=False, max_agents_by_species={'A': 100, 'B': 70, 'C': 140})

        self.sim1c = AbStochKin()
        self.sim1c.add_process_from_str("A + A -> 2C", k=0.005)
        self.sim1c.simulate(p0={'A': 70, 'C': 0}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 70, 'C': 70})

        self.sim1d = AbStochKin()
        self.sim1d.add_process_from_str("A + B -> C + D", k=0.01)
        self.sim1d.simulate(p0={'A': 30, 'B': 50, 'C': 0, 'D': 0}, t_max=10, dt=0.01, n=100,
                            show_plots=False,
                            max_agents_by_species={'A': 30, 'B': 50, 'C': 50, 'D': 50})

        self.sim1e = AbStochKin()
        self.sim1e.add_process_from_str("2A -> C + D", k=0.01)
        self.sim1e.simulate(p0={'A': 70, 'C': 0, 'D': 0}, t_max=10, dt=0.01, n=100,
                            show_plots=False, max_agents_by_species={'A': 70, 'C': 70, 'D': 70})

    def test_simulations_1(self):
        self.assertEqual(self.sim1a.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(unique(self.sim1a.sims[0].k_vals['A + B -> C, k = 0.03'])),
                            {0.03})
        self.assertGreaterEqual(self.sim1a.sims[0].results['A']['R^2'], 0.998)
        self.assertGreaterEqual(self.sim1a.sims[0].results['B']['R^2'], 0.998)
        self.assertGreaterEqual(self.sim1a.sims[0].results['C']['R^2'], 0.998)

        self.assertEqual(self.sim1b.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(unique(self.sim1b.sims[0].k_vals['A + B -> 2 C, k = 0.003'])),
                            {0.003})
        self.assertGreaterEqual(self.sim1b.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1b.sims[0].results['B']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1b.sims[0].results['C']['R^2'], 0.999)

        self.assertEqual(self.sim1c.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(unique(self.sim1c.sims[0].k_vals['2 A -> 2 C, k = 0.005'])),
                            {0.0, 0.005})
        self.assertGreaterEqual(self.sim1c.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1c.sims[0].results['C']['R^2'], 0.999)

        self.assertEqual(self.sim1d.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(unique(self.sim1d.sims[0].k_vals['A + B -> C + D, k = 0.01'])),
                            {0.01})
        self.assertGreaterEqual(self.sim1d.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1d.sims[0].results['B']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1d.sims[0].results['C']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1d.sims[0].results['D']['R^2'], 0.999)

        self.assertEqual(self.sim1e.sims[0]._het_processes_num, 0)
        self.assertSetEqual(set(unique(self.sim1e.sims[0].k_vals['2 A -> C + D, k = 0.01'])),
                            {0.0, 0.01})
        self.assertEqual(sum(diagonal(self.sim1e.sims[0].k_vals['2 A -> C + D, k = 0.01'])), 0)
        self.assertGreaterEqual(self.sim1e.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1e.sims[0].results['C']['R^2'], 0.998)
        self.assertGreaterEqual(self.sim1e.sims[0].results['D']['R^2'], 0.998)


if __name__ == '__main__':
    unittest.main()
