#  Copyright (c) 2024-2025, Alex Plakantonakis.
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

import os
import unittest

from abstochkin.base import AbStochKin


class TestAbStochKin(unittest.TestCase):
    def setUp(self):
        # Test importing processes from file 1
        self.sim1 = AbStochKin()
        self.sim1.add_processes_from_file(
            os.path.join(os.path.dirname(__file__), "processes_test_1.txt")
        )
        self.sim1.add_process({'G0_3': 1}, {'None': 0}, 1.023)
        self.sim1.add_process_from_str(" -> X", 1)
        self.sim1.simulate(p0={'C': 10, 'CaM_4Ca': 15, 'W_2': 40, 'Pi': 6, 'H2O': 100,
                               'W_1': 10, 'G0_3': 3, 'D': 15, 'W': 36, 'ATP': 22,
                               'ADP': 13, 'X': 57, 'Ca': 60, 'Y': 46, 'AMP': 14,
                               'CaM': 5, 'PPi': 0}, t_max=10, dt=0.1, n=100, solve_odes=False,
                           run=False)

        # Import two processes, each from a str
        self.sim2 = AbStochKin()
        self.sim2.add_process_from_str("2A -> B", 0.3)
        self.sim2.add_process_from_str("B -> ", 0.1)
        self.sim2.simulate(p0={'A': 100, 'B': 0}, t_max=10, dt=0.01, n=100, solve_odes=True,
                           run=False)

        # Test importing processes from file 2
        self.sim3 = AbStochKin()
        self.sim3.add_processes_from_file(
            os.path.join(os.path.dirname(__file__), "processes_test_2.txt")
        )

        # Test adding processes where the system is in a compartment with a specified volume
        self.sim4 = AbStochKin(volume=1.5e-15)  # Approximate volume of an E. coli cell
        self.sim4.add_process_from_str('2A <-> X', k=0.01, k_rev=0.05)
        self.sim4.add_process({'': 0}, {'C': 1},
                              k=0.001,
                              regulating_species='E', alpha=2.5, K50=10, nH=2)

    def test_add_processes(self):
        self.assertEqual(len(self.sim1.sims[0].all_species), 17)
        self.assertSetEqual(self.sim1.sims[0].all_species,
                            {'C', 'CaM_4Ca', 'W_2', 'Pi', 'H2O', 'W_1', 'G0_3', 'D',
                             'W', 'ATP', 'ADP', 'X', 'Ca', 'Y', 'AMP', 'CaM', 'PPi'})
        self.assertEqual(len(self.sim1.processes), 7)

        self.assertEqual(len(self.sim1.sims[0]._procs_by_reactant), 9)
        self.assertNotIn('', self.sim1.sims[0]._procs_by_reactant)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_reactant['W_1']), 1)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_reactant['CaM']), 1)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_reactant['H2O']), 2)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_reactant['ATP']), 2)

        self.assertEqual(len(self.sim1.sims[0]._procs_by_product), 8)
        self.assertNotIn('', self.sim1.sims[0]._procs_by_product)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_product['Y']), 1)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_product['Pi']), 1)
        self.assertEqual(len(self.sim1.sims[0]._procs_by_product['Pi']), 1)

        self.assertEqual(len(self.sim3.processes), 5)

        self.assertEqual(self.sim3.processes[0].k, 0.01)
        self.assertEqual(self.sim3.processes[0].k_rev, 0.05)

        self.assertEqual(self.sim3.processes[1].k, [0.01, 0.02])
        self.assertEqual(self.sim3.processes[1].k_rev, (0.05, 0.01))

        self.assertEqual(self.sim3.processes[2].k, 0.4)
        self.assertEqual(self.sim3.processes[2].catalyst, 'E')
        self.assertEqual(self.sim3.processes[2].Km, 10)

        self.assertEqual(self.sim3.processes[3].k, (0.55, 0.1))
        self.assertEqual(self.sim3.processes[3].regulating_species, 'R')
        self.assertEqual(self.sim3.processes[3].alpha, 2.5)
        self.assertEqual(self.sim3.processes[3].K50, [30, 20])
        self.assertEqual(self.sim3.processes[3].nH, 3)

        self.assertEqual(self.sim3.processes[4].k, (0.55, 0.1))
        self.assertEqual(self.sim3.processes[4].regulating_species, ['A', 'C'])
        self.assertEqual(self.sim3.processes[4].alpha, [2.5, 1])
        self.assertEqual(self.sim3.processes[4].K50, [(30, 5), [20, 10]])
        self.assertEqual(self.sim3.processes[4].nH, [3, 2])

        # Make sure k has been converted to its microscopic value
        self.assertAlmostEqual(self.sim4.processes[0].k, 1.11e-11, places=2)
        self.assertEqual(self.sim4.processes[0].k_rev, 0.05)
        self.assertEqual(self.sim4.processes[1].K50, 9033211140.0)

    def test_remove_processes(self):
        self.sim1.del_process({'C': 1, 'D': 1}, {'Y': 1}, k=0.01)
        self.assertEqual(len(self.sim1.processes), 6)

        self.sim3.del_process_from_str("D<->F", k=[0.01, 0.02], k_rev=(0.05, 0.01))
        self.assertEqual(len(self.sim3.processes), 4)

        self.sim4.del_process({'': 0}, {'C': 1},
                              k=0.001,
                              regulating_species='E', alpha=2.5, K50=10, nH=2)
        self.assertEqual(len(self.sim4.processes), 1)


if __name__ == '__main__':
    unittest.main()
