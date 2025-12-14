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

from abstochkin.base import AbStochKin
import unittest


class TestRegulatedFirstOrderSimulation1(unittest.TestCase):
    def setUp(self):
        # Complete Repression by B: alpha=0
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A -> B", k=0.3,
                                        regulating_species='B', alpha=0, K50=10, nH=3)
        self.sim1a.simulate(p0={'A': 50, 'B': 0}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 70, 'B': 70})

        # No regulation: alpha=1
        self.sim1b = AbStochKin()
        self.sim1b.add_process_from_str("A -> B", k=0.3,
                                        regulating_species='B', alpha=1, K50=10, nH=3)
        self.sim1b.simulate(p0={'A': 40, 'B': 10}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 60, 'B': 60})

        # Activation by B: alpha=2
        self.sim1c = AbStochKin()
        self.sim1c.add_process_from_str("A -> B", k=0.3,
                                        regulating_species='B', alpha=2, K50=10, nH=3)
        self.sim1c.simulate(p0={'A': 40, 'B': 10}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 60, 'B': 60})

        # Activation by A: alpha=1.5
        self.sim1d = AbStochKin()
        self.sim1d.add_process_from_str("A -> ", k=0.3,
                                        regulating_species='A', alpha=1.5, K50=10, nH=3)
        self.sim1d.simulate(p0={'A': 40}, t_max=10, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 40})

        # Two regulators: A, C: # Repression by A, activation by C
        self.sim1e = AbStochKin()
        self.sim1e.add_process_from_str("A -> C", k=0.3,
                                        regulating_species='A, C', alpha=[0, 2],
                                        K50=[15, 5], nH=[1, 2])
        self.sim1e.simulate(p0={'A': 50, 'C': 0}, t_max=15, dt=0.01, n=100, show_plots=False,
                            max_agents_by_species={'A': 50, 'C': 50})

    def test_simulations_1(self):
        """ Test simulations of process with homogeneous populations. """
        self.assertEqual(self.sim1a.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1a.sims[0].results['A']['R^2'], 0.99)
        self.assertGreaterEqual(self.sim1a.sims[0].results['B']['R^2'], 0.99)

        self.assertEqual(self.sim1b.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1b.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1b.sims[0].results['B']['R^2'], 0.999)

        self.assertEqual(self.sim1c.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1c.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1c.sims[0].results['B']['R^2'], 0.999)

        self.assertEqual(self.sim1d.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1d.sims[0].results['A']['R^2'], 0.999)

        self.assertEqual(self.sim1e.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1e.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1e.sims[0].results['C']['R^2'], 0.999)


if __name__ == '__main__':
    unittest.main()
