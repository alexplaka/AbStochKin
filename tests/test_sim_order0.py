"""
Test running a simulation of a 0th order process:  -> A
Test running a simulation of a regulated 0th order process:  -> A
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

from numpy import sqrt

from abstochkin.base import AbStochKin
from abstochkin.utils import r_squared


class TestZerothOrderSimulation(unittest.TestCase):
    def setUp(self):
        self.sim = AbStochKin()
        self.sim.add_process_from_str(' -> A', 1)
        self.sim.simulate(p0={'A': 0}, t_max=10, dt=0.01, n=1000, show_plots=False)

    def test_simulation(self):
        self.assertEqual(self.sim.sims[0]._het_processes_num, 0)
        self.assertDictEqual(self.sim.sims[0].k_vals, {})

        # Make sure mean simulation trajectory agrees with DE/CME prediction
        self.assertGreaterEqual(self.sim.sims[0].results['A']['R^2'], 0.999)

        """ Make sure the standard deviation of the ensemble of simulation 
        trajectories agrees with CME prediction. """
        cme_std = sqrt(self.sim.sims[0].de_calcs.odes_sol.sol(self.sim.sims[0].time))
        r_sq_std = r_squared(self.sim.sims[0].results['A']['N_std'], cme_std)
        self.assertGreaterEqual(r_sq_std, 0.98)

        """ Make sure the coefficient of variation of the ensemble of simulation 
        trajectories agrees with the prediction for a Poisson process. """
        r_sq_eta = r_squared(self.sim.sims[0].results['A']['eta'],
                             self.sim.sims[0].results['A']['eta_p'])
        self.assertGreaterEqual(r_sq_eta, 0.999)


class TestRegulatedZerothOrderSimulation(unittest.TestCase):
    def setUp(self):
        # Activation by A
        self.sim1 = AbStochKin()
        self.sim1.add_process_from_str(' -> A', k=0.5,
                                       regulating_species='A', alpha=2, K50=10, nH=3)
        self.sim1.simulate(p0={'A': 0}, t_max=20, dt=0.01, n=500, show_plots=False)

        # No regulation
        self.sim2 = AbStochKin()
        self.sim2.add_process_from_str(' -> A', k=0.5,
                                       regulating_species='A', alpha=1, K50=10, nH=3)
        self.sim2.simulate(p0={'A': 0}, t_max=20, dt=0.01, n=500, show_plots=False)

        # Repression by A
        self.sim3 = AbStochKin()
        self.sim3.add_process_from_str(' -> A', k=0.5,
                                       regulating_species='A', alpha=0, K50=10, nH=2)
        self.sim3.simulate(p0={'A': 0}, t_max=20, dt=0.01, n=500, show_plots=False)

        # Two regulators: Activation by A, repression by B
        self.sim4 = AbStochKin()
        self.sim4.add_process_from_str("A -> B", k=0.25)
        self.sim4.add_process_from_str(" -> C", k=0.5,
                                       regulating_species='A, B', alpha=[2, 0.5],
                                       K50=[10, 5], nH=[2, 1])
        self.sim4.simulate(p0={'A': 50, 'B': 0, 'C': 0}, t_max=10, dt=0.01, n=150,
                           max_agents_by_species={'A': 50, 'B': 50, 'C': 50},
                           show_plots=False)

    def test_simulations(self):
        self.assertEqual(self.sim1.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1.sims[0].results['A']['R^2'], 0.999)

        self.assertEqual(self.sim2.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim2.sims[0].results['A']['R^2'], 0.999)

        self.assertEqual(self.sim3.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim3.sims[0].results['A']['R^2'], 0.999)

        self.assertEqual(self.sim4.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim4.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim4.sims[0].results['B']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim4.sims[0].results['C']['R^2'], 0.989)


if __name__ == '__main__':
    unittest.main()
