""" Test the `DEcalcs` class (used for deterministic calculations). """

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

import os
import unittest

from sympy import symbols, Float
import numpy as np

from abstochkin.base import AbStochKin
from abstochkin.de_calcs import DEcalcs
from abstochkin.process import Process


class TestDEcalcs1(unittest.TestCase):
    def setUp(self):
        self.A, self.B, self.C, self.D, = symbols('A, B, C, D')

        self.sim1a = AbStochKin()  # 0th order
        self.sim1a.add_process_from_str(" -> B", k=0.5)
        self.sim1a.simulate(p0={'B': 0}, t_max=10, run=False)

        self.sim1b = AbStochKin()  # 1st order conversion
        self.sim1b.add_process_from_str("A -> B", k=0.1)
        self.sim1b.simulate(p0={'A': 100, 'B': 0}, t_max=10, run=False)

        self.sim1c = AbStochKin()  # 1st order degradation
        self.sim1c.add_process_from_str("A -> ", k=0.25)
        self.sim1c.simulate(p0={'A': 100}, t_max=10, run=False)

        self.sim1d = AbStochKin()  # birth-death process
        self.sim1d.add_process_from_str(" -> A", k=0.5)
        self.sim1d.add_process_from_str("A -> ", k=0.1)
        self.sim1d.simulate(p0={'A': 100}, t_max=10, run=False)

        self.sim2a = AbStochKin()  # Homologous 2nd order conversion - dimerization
        self.sim2a.add_process_from_str("2A -> C", k=0.01)
        self.sim2a.simulate(p0={'A': 100, 'C': 0}, t_max=10, run=False)

        self.sim2b = AbStochKin()  # Heterologous 2nd order conversion - dimerization
        self.sim2b.add_process_from_str("A + B -> C", k=0.01)
        self.sim2b.simulate(p0={'A': 100, 'B': 100, 'C': 0}, t_max=10, run=False)

        self.sim2c = AbStochKin()  # Homologous 2nd order conversion
        self.sim2c.add_process_from_str("2A -> C + D", k=0.01)
        self.sim2c.simulate(p0={'A': 100, 'C': 0, 'D': 0}, t_max=10, run=False)

        self.sim2d = AbStochKin()  # Heterologous 2nd order conversion
        self.sim2d.add_process_from_str("A + B -> C + D", k=0.01)
        self.sim2d.simulate(p0={'A': 100, 'B': 100, 'C': 0, 'D': 0}, t_max=10, run=False)

    def test_setup_ODEs(self):
        self.assertDictEqual(self.sim1a.sims[0].de_calcs.odes, {'B': 0.5})
        self.assertDictEqual(self.sim1b.sims[0].de_calcs.odes,
                             {'A': -0.1 * self.A,
                              'B': 0.1 * self.A})
        self.assertDictEqual(self.sim1c.sims[0].de_calcs.odes,
                             {'A': -0.25 * self.A})
        self.assertDictEqual(self.sim1d.sims[0].de_calcs.odes,
                             {'A': 0.5 - 0.1 * self.A})

        self.assertDictEqual(self.sim2a.sims[0].de_calcs.odes,
                             {'A': -0.02 * self.A * (self.A - 1),
                              'C': 0.01 * self.A * (self.A - 1)})
        self.assertDictEqual(self.sim2b.sims[0].de_calcs.odes,
                             {'A': -0.01 * self.A * self.B,
                              'B': -0.01 * self.A * self.B,
                              'C': 0.01 * self.A * self.B})
        self.assertDictEqual(self.sim2c.sims[0].de_calcs.odes,
                             {'A': -0.02 * self.A * (self.A - 1),
                              'C': 0.01 * self.A * (self.A - 1),
                              'D': 0.01 * self.A * (self.A - 1)})
        self.assertDictEqual(self.sim2d.sims[0].de_calcs.odes,
                             {'A': -0.01 * self.A * self.B,
                              'B': -0.01 * self.A * self.B,
                              'C': 0.01 * self.A * self.B,
                              'D': 0.01 * self.A * self.B})

    def test_fixed_pts(self):
        self.assertEqual(True, True)


class TestDEcalcs2(unittest.TestCase):
    def setUp(self):
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

        self.sim2 = AbStochKin()
        self.sim2.add_process_from_str("2A -> B", 0.3)
        self.sim2.add_process_from_str("B -> ", 0.1)
        self.sim2.simulate(p0={'A': 100, 'B': 0}, t_max=10, dt=0.01, n=100, solve_odes=True,
                           run=False)

        # self.maxDiff = None

    def test_odes(self):
        # First, create sympy symbols for `sim1` object:
        X, Pi, PPi, CaM_4Ca, W_1, W_2, H2O, Y, C, D, CaM, W, G0_3, ATP, AMP, ADP, Ca = symbols(
            ['X', 'Pi', 'PPi', 'CaM_4Ca', 'W_1', 'W_2', 'H2O', 'Y', 'C', 'D', 'CaM', 'W', 'G0_3',
             'ATP', 'AMP', 'ADP', 'Ca'])
        # The following is the expected species-specific set of ODEs for `sim1`:
        odes1 = {'Y': 0.01 * C * D,
                 'W': 0.413 * W_1 * W_2,
                 'Pi': 0.3435 * ATP * H2O,
                 'ADP': 0.3435 * ATP * H2O,
                 'AMP': 0.4562 * ATP * H2O,
                 'PPi': 0.4562 * ATP * H2O,
                 'X': Float(1),
                 'D': -0.01 * C * D,
                 'C': -0.01 * C * D,
                 'W_2': -0.413 * W_1 * W_2,
                 'W_1': -0.413 * W_1 * W_2,
                 'CaM_4Ca': 0.15 * Ca * (Ca - 1) * (Ca - 2) * (Ca - 3) * CaM,
                 'CaM': -0.15 * Ca * (Ca - 1) * (Ca - 2) * (Ca - 3) * CaM,
                 'Ca': -0.6 * Ca * (Ca - 1) * (Ca - 2) * (Ca - 3) * CaM,
                 'H2O': -(0.3435 + 0.4562) * ATP * H2O,
                 'ATP': -(0.3435 + 0.4562) * ATP * H2O,
                 'G0_3': -1.023 * G0_3}
        self.assertDictEqual(self.sim1.sims[0].de_calcs.odes, odes1)

        """ Solution times out for this system of ODEs and initial conditions.
        Since this testing scenario was devised solely for testing the internal 
        data structures, we do not test whether the solution was reached. """
        # self.assertTrue(self.sim1.sims[0].odes_sol.success)
        # self.assertEqual(len(self.sim1.sims[0].odes_sol.y), len(self.sim1.all_species))

        A, B = symbols(['A', 'B'])  # Create sympy symbols for `sim2` object
        odes2 = {'A': -0.6 * A * (A - 1),
                 'B': 0.3 * A * (A - 1) - 0.1 * B}  # expected ODEs for `sim2`
        self.assertDictEqual(self.sim2.sims[0].de_calcs.odes, odes2)
        self.assertTrue(self.sim2.sims[0].de_calcs.odes_sol.success)
        self.assertEqual(len(self.sim2.sims[0].de_calcs.odes_sol.y),
                         len(self.sim2.sims[0].all_species))


class TestDEcalcsMM(unittest.TestCase):
    """ Test that the ODEs for a Michaelis-Menten Process are correctly set up. """

    def setUp(self):
        self.sim = AbStochKin()
        self.sim.add_process_from_str("A -> B", 0.3, catalyst='E', Km=5)
        self.sim.simulate(p0={'A': 40, 'B': 0, 'E': 10}, t_max=20, solve_odes=True, run=False)

    def test_odes_mm(self):
        A, B, E = symbols('A, B, E')
        self.assertEqual(self.sim.sims[0].de_calcs.odes['A'], -0.3 * E * A / (A + 5.0))
        self.assertEqual(self.sim.sims[0].de_calcs.odes['B'], 0.3 * E * A / (A + 5.0))
        self.assertEqual(self.sim.sims[0].de_calcs.odes['E'], 0)


class TestDEcalcsReg(unittest.TestCase):
    """ Test that the ODEs for a Regulated Process are correctly set up. """

    def setUp(self):
        """
        Set up two regulated processes:
        1. A process regulated by one species.
        2. A process regulated by two species.
        """
        self.sim1 = AbStochKin()
        self.sim1.add_process_from_str('A->B', k=0.3,
                                       regulating_species='B', alpha=2, K50=10, nH=3)
        self.sim1.simulate(p0={'A': 40, 'B': 10}, t_max=10, solve_odes=True, run=False)

        self.sim2 = AbStochKin()
        self.sim2.add_process_from_str('A->B', k=0.3,
                                       regulating_species='A, B', alpha=[2, 0],
                                       K50=[10, (20, 5)], nH=[3, 2])
        self.sim2.simulate(p0={'A': 40, 'B': 10}, t_max=10, solve_odes=True, run=False)

    def test_odes_reg(self):
        A, B = symbols('A, B')

        expr1 = 0.3 * A * (1 + 2 * (0.1 * B) ** 3) / (1 + (0.1 * B) ** 3)
        self.assertEqual(self.sim1.sims[0].de_calcs.odes['A'], -expr1)
        self.assertEqual(self.sim1.sims[0].de_calcs.odes['B'], expr1)

        reg_term_2a = (1 + 2 * (0.1 * A) ** 3) / (1 + (0.1 * A) ** 3)
        reg_term_2b = 1 / (1 + (B / 20) ** 2)
        expr2 = 0.3 * A * reg_term_2a * reg_term_2b
        self.assertEqual(self.sim2.sims[0].de_calcs.odes['A'], -expr2)
        self.assertEqual(self.sim2.sims[0].de_calcs.odes['B'], expr2)


class TestDEcalcsRegMM(unittest.TestCase):
    """ Test that the ODEs for a Regulated Michaelis-Menten Process
    are correctly set up. """

    def setUp(self):
        self.sim1 = AbStochKin()
        self.sim1.add_process_from_str('A->B', 0.3,
                                       regulating_species='A, B',
                                       alpha=[2, 0.5],
                                       K50=[20, 10],
                                       nH=[1, 3],
                                       catalyst='C',
                                       Km=5)
        self.sim1.simulate(p0={'A': 50, 'B': 0, 'C': 10}, t_max=10, dt=0.01, n=100, run=False,
                           max_agents_by_species={'A': 50, 'B': 0, 'C': 10})

    def test_odes_reg_mm(self):
        A, B, C = symbols('A, B, C')

        reg_term_1a = (1 + 2 * 0.05 * A) / (1 + 0.05 * A)
        reg_term_1b = (1 + 0.5 * (0.1 * B) ** 3) / (1 + (0.1 * B) ** 3)
        reg_term_1c = C / (A + 5.0)
        expr1 = 0.3 * A * reg_term_1a * reg_term_1b * reg_term_1c
        self.assertEqual(self.sim1.sims[0].de_calcs.odes['A'], -expr1)
        self.assertEqual(self.sim1.sims[0].de_calcs.odes['B'], expr1)
        self.assertEqual(self.sim1.sims[0].de_calcs.odes['C'], 0)


class TestGetFixedPtsNumerically(unittest.TestCase):
    """ Test that the fixed points are correctly found numerically. """

    def setUp(self):
        self.obj0 = DEcalcs(p0={'A': 20},
                            t_min=0,
                            t_max=5,
                            processes=[Process.from_string("A -> ", k=0.1)],  # dummy process
                            ode_method="RK45",
                            time_unit='sec')
        # self.obj0._odes_lambdified = lambda x: np.array([x[0] - 1, x[1] - 2])
        self.obj0.get_fixed_pts_numerically()

    def test_success(self):
        self.assertTrue(self.obj0.fixed_pts_num_sol_info.success)

    def test_fixed_pts_num_sol(self):
        np.testing.assert_equal(self.obj0.fixed_pts_num_sol['A'], 0)

    def test_root_method(self):
        self.obj0.get_fixed_pts_numerically(root_method='lm')
        self.assertEqual(self.obj0.fixed_pts_num_sol_info.method, 'lm')

    def test_failure(self):
        pass


if __name__ == '__main__':
    unittest.main()
