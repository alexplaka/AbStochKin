"""
Test running a simulation of two birth-death processes:
    - -> A ->
    - -> A -> B
In both cases, the 1st process of A is assumed to be homogeneous.
This way, the statistics of the simulated trajectories can be
compared to deterministic (ODE and CME) predictions.
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

from abstochkin import AbStochKin
from abstochkin.utils import r_squared


class BirthDeathProcess1(unittest.TestCase):
    def setUp(self):
        self.sim = AbStochKin()
        self.sim.add_process_from_str(" -> A", k=2)
        self.sim.add_process_from_str("A -> ", k=0.4)

    def test_simulation(self):
        self.sim.simulate(p0={'A': 0}, t_max=10, dt=0.01, n=250, show_plots=False,
                          max_agents_by_species={'A': 20})

        self.assertEqual(self.sim.sims[0]._het_processes_num, 0)

        # Make sure mean simulation trajectory agrees with DE/CME prediction
        self.assertGreaterEqual(self.sim.sims[0].results['A']['R^2'], 0.99)

        """ Make sure the standard deviation of the ensemble of simulation 
        trajectories agrees with CME prediction = sqrt(<A>). """
        cme_std = sqrt(self.sim.sims[0].de_calcs.odes_sol.sol(self.sim.sims[0].time))
        r_sq_std = r_squared(self.sim.sims[0].results['A']['N_std'], cme_std)
        self.assertGreaterEqual(r_sq_std, 0.85)
        """ This is a very noisy process. Would need thousands of repetitions 
        of the simulation to get `r_sq_std` to be greater than `0.99`. This is 
        why we set the threshold for the test to pass to be only `0.85`. """

        """ However, the coefficient of variation gives the standard deviation 
        adjusted by the mean population trajectory, so it's a more 
        reliable indicator of the level of noise.  
        Make sure the coefficient of variation of the ensemble of simulation 
        trajectories agrees with the prediction for a Poisson process. """
        r_sq_eta = r_squared(self.sim.sims[0].results['A']['eta'],
                             self.sim.sims[0].results['A']['eta_p'])
        self.assertGreaterEqual(r_sq_eta, 0.99)


class BirthDeathProcess2(unittest.TestCase):
    def setUp(self):
        self.sim = AbStochKin()
        self.sim.add_process_from_str(" -> A", k=1)
        self.sim.add_process_from_str("A -> B", k=0.15)

    def test_simulation(self):
        self.sim.simulate(p0={'A': 0, 'B': 0}, t_max=10, dt=0.01, n=250, show_plots=False)
        self.assertEqual(self.sim.sims[0]._het_processes_num, 0)

        # Make sure mean simulation trajectory agrees with DE/CME prediction
        self.assertGreaterEqual(self.sim.sims[0].results['A']['R^2'], 0.989)
        self.assertGreaterEqual(self.sim.sims[0].results['B']['R^2'], 0.999)

        for sp in ['A', 'B']:
            self.assertGreaterEqual(
                r_squared(self.sim.sims[0].results[sp]['eta'],
                          self.sim.sims[0].results[sp]['eta_p']),
                0.99
            )


if __name__ == '__main__':
    unittest.main()
