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


class TestRegulatedMichaelisMentenFirstOrderSimulation(unittest.TestCase):
    def setUp(self):
        self.sim1 = AbStochKin()
        self.sim1.add_process_from_str('A->B', 0.3,
                                       regulating_species='A', alpha=2,
                                       K50=10, nH=1,
                                       catalyst='C', Km=5)
        self.sim1.simulate(p0={'A': 50, 'B': 0, 'C': 10}, t_max=20, dt=0.01, n=100,
                           show_plots=False, multithreading=False,
                           max_agents_by_species={'A': 50, 'B': 50, 'C': 50})

    def test_simulations(self):
        self.assertEqual(self.sim1.sims[0]._het_processes_num, 0)
        self.assertGreaterEqual(self.sim1.sims[0].results['A']['R^2'], 0.999)
        self.assertGreaterEqual(self.sim1.sims[0].results['B']['R^2'], 0.999)


if __name__ == '__main__':
    unittest.main()
