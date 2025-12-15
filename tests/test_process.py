""" Test the Process class and its subclasses. """

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

from abstochkin.process import NullSpeciesNameError
from abstochkin.process import Process, ReversibleProcess, MichaelisMentenProcess, \
    RegulatedProcess, RegulatedMichaelisMentenProcess


class TestProcess(unittest.TestCase):
    def setUp(self):
        """ Set up some processes for testing. """
        # Some processes from a dictionary:
        self.proc1a = Process({'A': 2}, {'B': 1}, 0.01)
        self.proc1b = Process({'A': 1, 'M_2': 1}, {'W': 2}, 0.1)
        self.proc1c = Process({'': 0}, {'A_d': 1}, 0.157)
        self.proc1d = Process({'G0-3': 1}, {'None': 0}, 1.023)

        # Some processes from a string:
        self.proc2a = Process.from_string('X + 3 W --> 2X', 0.001)
        self.proc2b = Process.from_string(' -> 2 X', 0.349)
        self.proc2c = Process.from_string('A ,  ', 0.45, sep=',')
        self.proc2d = Process.from_string('X_protein ->  None  ', 1)
        self.proc2e = Process.from_string('A + A ->  C  ', 0.5)
        self.proc2f = Process.from_string('2A ->  C  ', 0.5)
        self.proc2g = Process.from_string('A -> 2 C + C  ', 0.05)

        # Some other processes:
        self.proc3a = Process.from_string('-> A', 1)
        self.proc3b = Process.from_string('2Yto3 z', 1e-2, sep='to')
        self.proc3c = Process.from_string('2A 4G-->5 F', 2E-4)

    def test_reactants(self):
        self.assertDictEqual(self.proc1a.reactants, {'A': 2})
        self.assertDictEqual(self.proc1b.reactants, {'A': 1, 'M_2': 1})
        self.assertDictEqual(self.proc1c.reactants, {'': 0})
        self.assertDictEqual(self.proc1d.reactants, {'G0-3': 1})

        self.assertDictEqual(self.proc2a.reactants, {'X': 1, 'W': 3})
        self.assertDictEqual(self.proc2b.reactants, {'': 0})
        self.assertDictEqual(self.proc2c.reactants, {'A': 1})
        self.assertDictEqual(self.proc2d.reactants, {'X_protein': 1})
        self.assertDictEqual(self.proc2e.reactants, {'A': 2})
        self.assertDictEqual(self.proc2f.reactants, {'A': 2})
        self.assertDictEqual(self.proc2g.reactants, {'A': 1})

        self.assertDictEqual(self.proc3a.reactants, {'': 0})
        self.assertDictEqual(self.proc3b.reactants, {'Y': 2})
        self.assertDictEqual(self.proc3c.reactants, {'A 4G': 2})
        self.assertRaises(NullSpeciesNameError, Process.from_string, '2->', 0.56)

        self.assertEqual(self.proc1a.reacts_, ['A', 'A'])
        self.assertEqual(self.proc1b.reacts_, ['A', 'M_2'])
        self.assertEqual(self.proc1c.reacts_, [])
        self.assertEqual(self.proc1d.reacts_, ['G0-3'])

        self.assertEqual(self.proc2a.reacts_, ['X', 'W', 'W', 'W'])
        self.assertEqual(self.proc2b.reacts_, [])
        self.assertEqual(self.proc2c.reacts_, ['A'])
        self.assertEqual(self.proc2d.reacts_, ['X_protein'])
        self.assertEqual(self.proc2e.reacts_, ['A', 'A'])
        self.assertEqual(self.proc2f.reacts_, ['A', 'A'])
        self.assertEqual(self.proc2g.reacts_, ['A'])

        self.assertEqual(self.proc3a.reacts_, [])
        self.assertEqual(self.proc3b.reacts_, ['Y', 'Y'])
        self.assertEqual(self.proc3c.reacts_, ['A 4G', 'A 4G'])

    def test_products(self):
        self.assertDictEqual(self.proc1a.products, {'B': 1})
        self.assertDictEqual(self.proc1b.products, {'W': 2})
        self.assertDictEqual(self.proc1c.products, {'A_d': 1})
        self.assertDictEqual(self.proc1d.products, {'': 0})

        self.assertDictEqual(self.proc2a.products, {'X': 2})
        self.assertDictEqual(self.proc2b.products, {'X': 2})
        self.assertDictEqual(self.proc2c.products, {'': 0})
        self.assertDictEqual(self.proc2d.products, {'': 0})
        self.assertDictEqual(self.proc2e.products, {'C': 1})
        self.assertDictEqual(self.proc2f.products, {'C': 1})
        self.assertDictEqual(self.proc2g.products, {'C': 3})

        self.assertDictEqual(self.proc3a.products, {'A': 1})
        self.assertDictEqual(self.proc3b.products, {'z': 3})
        self.assertDictEqual(self.proc3c.products, {'F': 5})

        # output = io.StringIO()
        # handler = logging.StreamHandler(output)
        # logger = logging.getLogger()
        # logger.addHandler(handler)
        # with redirect_stdout(output):
        #     Process.from_string('2A 4G-->-5 F', 2E-4)
        # self.assertEqual(output.getvalue(), "Assertion failed: Coefficient cannot be negative: -5, species F.\n")

        self.assertRaises(AssertionError, Process.from_string, '2A 4G-->-5 F', 2E-4)

        self.assertRaises(NullSpeciesNameError, Process.from_string, '2A + X -> 3', 0.11)

        self.assertEqual(self.proc1a.prods_, ['B'])
        self.assertEqual(self.proc1b.prods_, ['W', 'W'])
        self.assertEqual(self.proc1c.prods_, ['A_d'])
        self.assertEqual(self.proc1d.prods_, [])

        self.assertEqual(self.proc2a.prods_, ['X', 'X'])
        self.assertEqual(self.proc2b.prods_, ['X', 'X'])
        self.assertEqual(self.proc2c.prods_, [])
        self.assertEqual(self.proc2d.prods_, [])
        self.assertEqual(self.proc2e.prods_, ['C'])
        self.assertEqual(self.proc2f.prods_, ['C'])
        self.assertEqual(self.proc2g.prods_, ['C', 'C', 'C'])

        self.assertEqual(self.proc3a.prods_, ['A'])
        self.assertEqual(self.proc3b.prods_, ['z', 'z', 'z'])
        self.assertEqual(self.proc3c.prods_, ['F', 'F', 'F', 'F', 'F'])

    def test_order(self):
        self.assertEqual(self.proc1a.order, 2)
        self.assertEqual(self.proc1b.order, 2)
        self.assertEqual(self.proc1c.order, 0)
        self.assertEqual(self.proc1d.order, 1)

        self.assertEqual(self.proc2a.order, 4)
        self.assertEqual(self.proc2b.order, 0)
        self.assertEqual(self.proc2c.order, 1)
        self.assertEqual(self.proc2d.order, 1)
        self.assertEqual(self.proc2e.order, 2)
        self.assertEqual(self.proc2f.order, 2)
        self.assertEqual(self.proc2e, self.proc2f)
        self.assertEqual(self.proc2g.order, 1)

        self.assertEqual(self.proc3a.order, 0)
        self.assertEqual(self.proc3b.order, 2)
        self.assertEqual(self.proc3c.order, 2)


class TestMichaelisMentenProcess(unittest.TestCase):
    def setUp(self):
        self.proc1a = MichaelisMentenProcess.from_string('A --> B', 0.001, catalyst='E', Km=15)

    def test_MM_processes(self):
        self.assertEqual(self.proc1a.catalyst, 'E')
        self.assertEqual(self.proc1a.Km, 15)
        self.assertIs(self.proc1a.is_heterogeneous, False)
        self.assertIs(self.proc1a.is_heterogeneous_Km, False)
        self.assertSetEqual(self.proc1a.species, {'A', 'B', 'E'})


class TestReversibleProcess(unittest.TestCase):
    def setUp(self):
        self.proc1a = ReversibleProcess.from_string('A <--> B', 0.1, k_rev=0.2)
        self.proc1b = ReversibleProcess.from_string('A <-> 2B',
                                                    k=[0.5, 0.25], k_rev=(0.3, 0.1))

    def test_reversible_processes(self):
        self.assertEqual(self.proc1a.k, 0.1)
        self.assertEqual(self.proc1a.k_rev, 0.2)
        self.assertEqual(self.proc1a.order, 1)
        self.assertEqual(self.proc1a.order_rev, 1)
        self.assertEqual(self.proc1a.is_heterogeneous, False)
        self.assertEqual(self.proc1a.is_heterogeneous_rev, False)

        self.assertListEqual(self.proc1b.k, [0.5, 0.25])
        self.assertTupleEqual(self.proc1b.k_rev, (0.3, 0.1))
        self.assertEqual(self.proc1b.order, 1)
        self.assertEqual(self.proc1b.order_rev, 2)
        self.assertEqual(self.proc1b.is_heterogeneous, True)
        self.assertEqual(self.proc1b.is_heterogeneous_rev, True)


class TestRegulatedProcess(unittest.TestCase):
    def setUp(self):
        self.proc1a = RegulatedProcess.from_string('A --> B', 0.15,
                                                   regulating_species='A', alpha=2, K50=5, nH=1)
        self.proc1b = RegulatedProcess.from_string('A + B -> C',
                                                   k=[0.5, 0.25], regulating_species='A, C ',
                                                   alpha=[2, 0], K50=[10, (5, 1)], nH=[3, 2])

    def test_regulated_processes(self):
        self.assertEqual(self.proc1a.k, 0.15)
        self.assertEqual(self.proc1a.is_heterogeneous, False)
        self.assertEqual(self.proc1a.order, 1)
        self.assertEqual(self.proc1a.regulating_species, 'A')
        self.assertEqual(self.proc1a.is_heterogeneous_K50, False)
        self.assertEqual(self.proc1a.alpha, 2)
        self.assertEqual(self.proc1a.K50, 5)
        self.assertEqual(self.proc1a.nH, 1)
        self.assertEqual(self.proc1a.regulation_type, 'activation')

        self.assertEqual(self.proc1b.k, [0.5, 0.25])
        self.assertEqual(self.proc1b.is_heterogeneous, True)
        self.assertEqual(self.proc1b.order, 2)
        self.assertEqual(self.proc1b.regulating_species, ['A', 'C'])
        self.assertEqual(self.proc1b.is_heterogeneous_K50, [False, True])
        self.assertEqual(self.proc1b.alpha, [2, 0])
        self.assertEqual(self.proc1b.K50, [10, (5, 1)])
        self.assertEqual(self.proc1b.nH, [3, 2])
        self.assertEqual(self.proc1b.regulation_type, ['activation', 'repression'])


class TestRegulatedMichaelisMentenProcess(unittest.TestCase):
    def setUp(self):
        self.proc1a = RegulatedMichaelisMentenProcess.from_string('A --> B', 0.05,
                                                                  regulating_species='A', alpha=10,
                                                                  K50=5, nH=1,
                                                                  catalyst='C', Km=20)
        self.proc1b = RegulatedMichaelisMentenProcess.from_string('X -> Y',
                                                                  k=[0.5, 0.25],
                                                                  regulating_species=' X , Y ',
                                                                  alpha=[2, 0],
                                                                  K50=[(5, 1), 10], nH=[3, 2],
                                                                  catalyst='E', Km=(10, 2))

    def test_regulated_mm_processes(self):
        self.assertEqual(self.proc1a.k, 0.05)
        self.assertEqual(self.proc1a.is_heterogeneous, False)
        self.assertEqual(self.proc1a.order, 1)
        self.assertEqual(self.proc1a.regulating_species, 'A')
        self.assertEqual(self.proc1a.is_heterogeneous_K50, False)
        self.assertEqual(self.proc1a.alpha, 10)
        self.assertEqual(self.proc1a.K50, 5)
        self.assertEqual(self.proc1a.nH, 1)
        self.assertEqual(self.proc1a.regulation_type, 'activation')
        self.assertEqual(self.proc1a.catalyst, 'C')
        self.assertEqual(self.proc1a.Km, 20)
        self.assertEqual(self.proc1a.is_heterogeneous_Km, False)
        self.assertEqual(self.proc1a.species, {'A', 'B', 'C'})

        self.assertEqual(self.proc1b.k, [0.5, 0.25])
        self.assertEqual(self.proc1b.is_heterogeneous, True)
        self.assertEqual(self.proc1b.order, 1)
        self.assertEqual(self.proc1b.regulating_species, ['X', 'Y'])
        self.assertEqual(self.proc1b.is_heterogeneous_K50, [True, False])
        self.assertEqual(self.proc1b.alpha, [2, 0])
        self.assertEqual(self.proc1b.K50, [(5, 1), 10])
        self.assertEqual(self.proc1b.nH, [3, 2])
        self.assertEqual(self.proc1b.regulation_type, ['activation', 'repression'])
        self.assertEqual(self.proc1b.catalyst, 'E')
        self.assertEqual(self.proc1b.Km, (10, 2))
        self.assertEqual(self.proc1b.is_heterogeneous_Km, True)
        self.assertEqual(self.proc1b.species, {'X', 'Y', 'E'})


class TestProcessEquality(unittest.TestCase):
    def setUp(self):
        self.a = Process.from_string('A->B', k=0.1)
        self.b = ReversibleProcess.from_string('A<->B', k=0.1, k_rev=0.2)
        self.c = MichaelisMentenProcess.from_string('A->B', k=0.1, catalyst='C', Km=10)
        self.d = RegulatedProcess.from_string('A->B', k=0.1, regulating_species='B',
                                              alpha=2, K50=10, nH=3)
        self.e = RegulatedProcess.from_string('A->B', k=0.1, regulating_species='A, B ',
                                              alpha=[2, 1], K50=[10, 5], nH=[3, 2])
        self.f = RegulatedMichaelisMentenProcess.from_string('A->B', k=0.1,
                                                             regulating_species='A, B ',
                                                             alpha=[2, 0], K50=[10, 5], nH=[3, 2],
                                                             catalyst='E', Km=12.5)

        self.procs_dict = dict()
        self.procs_dict[self.a] = 0
        self.procs_dict[self.b] = 1
        self.procs_dict[self.c] = 2
        self.procs_dict[self.d] = 3
        self.procs_dict[self.e] = 4
        self.procs_dict[self.f] = 5

    def test_procs_eq(self):
        # Check that they are not equal to each other.
        self.assertEqual(len({self.a, self.b, self.c, self.d, self.e}), 5)

        # Check that the hashes are not equal to each other.
        self.assertEqual(len({hash(self.a), hash(self.b), hash(self.c), hash(self.d), hash(self.e)}), 5)

        a_str = 'A -> B, k = 0.1'
        b_str = 'A <-> B, k = 0.1, k_rev = 0.2'
        c_str = 'A -> B, k = 0.1, catalyst = C, Km = 10'
        d_str = 'A -> B, k = 0.1, regulating_species = B, alpha = 2, K50 = 10, nH = 3'
        e_str = "A -> B, k = 0.1, regulating_species = ['A', 'B'], alpha = [2, 1], " \
                "K50 = [10, 5], nH = [3, 2]"
        f_str = "A -> B, k = 0.1, regulating_species = ['A', 'B'], alpha = [2, 0], " \
                "K50 = [10, 5], nH = [3, 2], catalyst = E, Km = 12.5"

        # Check short string representation of each process
        self.assertEqual(self.a._str, a_str)
        self.assertEqual(self.b._str, b_str)
        self.assertEqual(self.c._str, c_str)
        self.assertEqual(self.d._str, d_str)
        self.assertEqual(self.e._str, e_str)
        self.assertEqual(self.f._str, f_str)

        # Check equality of a process to a string
        self.assertEqual(self.a, a_str)
        self.assertEqual(self.b, b_str)
        self.assertEqual(self.c, c_str)
        self.assertEqual(self.d, d_str)
        self.assertEqual(self.e, e_str)
        self.assertEqual(self.f, f_str)

        # Check equality of a process to a string without whitespaces
        self.assertEqual(self.a, a_str.replace(' ', ''))
        self.assertEqual(self.b, b_str.replace(' ', ''))
        self.assertEqual(self.c, c_str.replace(' ', ''))
        self.assertEqual(self.d, d_str.replace(' ', ''))
        self.assertEqual(self.e, e_str.replace(' ', ''))
        self.assertEqual(self.f, f_str.replace(' ', ''))

        # Check equality of a process to itself
        self.assertEqual(self.a, self.a)
        self.assertEqual(self.b, self.b)
        self.assertEqual(self.c, self.c)
        self.assertEqual(self.d, self.d)
        self.assertEqual(self.e, self.e)
        self.assertEqual(self.f, self.f)

        # Check referring to a process as a dictionary key
        self.assertEqual(self.procs_dict[self.a], self.procs_dict[a_str])
        self.assertEqual(self.procs_dict[self.b], self.procs_dict[b_str])
        self.assertEqual(self.procs_dict[self.c], self.procs_dict[c_str])
        self.assertEqual(self.procs_dict[self.d], self.procs_dict[d_str])
        self.assertEqual(self.procs_dict[self.e], self.procs_dict[e_str])
        self.assertEqual(self.procs_dict[self.f], self.procs_dict[f_str])


class TestProcessInVolume(unittest.TestCase):
    def setUp(self):
        """ Set up some processes for testing. """
        # Some processes from a dictionary:
        self.proc1a = Process({'A': 2}, {'B': 1}, 0.01, volume=1e-9)
        self.proc1b = Process({'A': 1}, {'B': 1}, 0.01, volume=1e-9)
        self.proc1c = Process({'': 0}, {'A_d': 1}, 0.157, volume=1e-9)

        # Some processes from a string:
        self.proc2a = Process.from_string(' -> 2 X', 0.349, volume=1e-6)
        self.proc2b = Process.from_string('A + A ->  C  ', 0.5, volume=1e-6)
        self.proc2c = Process.from_string('A -> 2 C + C  ', (0.05, 0.01), volume=1e-6)

    def test_micro_k_vals(self):
        self.assertAlmostEqual(self.proc1a.k, 1.6605390671738466e-17)
        self.assertEqual(self.proc1b.k, 0.01)
        self.assertAlmostEqual(self.proc1c.k, 94547609932000.0)

        self.assertAlmostEqual(self.proc2a.k, 2.10172712524e+17)
        self.assertAlmostEqual(self.proc2b.k, 8.302695335869233e-19)
        self.assertTupleEqual(self.proc2c.k, (0.05, 0.01))


class TestReversibleProcessInVolume(unittest.TestCase):
    def setUp(self):
        self.proc1a = ReversibleProcess({'A': 1}, {'B': 1},
                                        k=0.1, k_rev=0.2,
                                        volume=1e-10)
        self.proc1b = ReversibleProcess({'A': 1}, {'B': 2},
                                        k=[0.5, 0.25, 0.35], k_rev=(0.3, 0.1),
                                        volume=1e-8)

        self.proc2a = ReversibleProcess.from_string('A <--> B',
                                                    k=0.2, k_rev=0.15,
                                                    volume=3e-9)
        self.proc2b = ReversibleProcess.from_string('A <-> 2B',
                                                    k=[0.11, 0.29], k_rev=[0.25, 0.05],
                                                    volume=2e-8)
        self.proc2c = ReversibleProcess.from_string('2A <-> B',
                                                    k=(0.5, 0.25), k_rev=0.1,
                                                    volume=2e-10)

    def test_micro_k_vals(self):
        self.assertEqual(self.proc1a.k, 0.1)
        self.assertEqual(self.proc1a.k_rev, 0.2)
        self.assertListEqual(self.proc1b.k, [0.5, 0.25, 0.35])
        self.assertTupleEqual(self.proc1b.k_rev, (4.98161720152154e-17, 1.660539067173847e-17))

        self.assertEqual(self.proc2a.k, 0.2)
        self.assertEqual(self.proc2a.k_rev, 0.15)
        self.assertListEqual(self.proc2b.k, [0.11, 0.29])
        self.assertListEqual(self.proc2b.k_rev, [2.0756738339673084e-17, 4.151347667934617e-18])

        self.assertTupleEqual(self.proc2c.k, (4.1513476679346165e-15, 2.0756738339673083e-15))
        self.assertEqual(self.proc2c.k_rev, 0.1)


class TestMichaelisMentenProcessInVolume(unittest.TestCase):
    def setUp(self):
        self.proc1a = MichaelisMentenProcess({'A': 1}, {'B': 1},
                                             k=0.03, catalyst='E', Km=1e-5,
                                             volume=1e-15)
        self.proc1b = MichaelisMentenProcess({'A': 1}, {'B': 1},
                                             k=(0.03, 0.01), catalyst='E',
                                             Km=[1e-5, 2e-5],
                                             volume=1e-15)

        self.proc1c = MichaelisMentenProcess.from_string("X -> Y",
                                                         k=[0.05, 0.04, 0.03],
                                                         catalyst='E',
                                                         Km=(2e-8, 0.5e-8),
                                                         volume=1e-16)

    def test_micro_k_vals(self):
        self.assertEqual(self.proc1a.k, 0.03)
        self.assertAlmostEqual(self.proc1a.Km, 6022.14076)
        self.assertEqual(self.proc1b.k, (0.03, 0.01))
        self.assertIsInstance(self.proc1b.Km, list)
        self.assertAlmostEqual(self.proc1b.Km[0], 6022.14076)
        self.assertAlmostEqual(self.proc1b.Km[1], 2 * 6022.14076)

        self.assertListEqual(self.proc1c.k, [0.05, 0.04, 0.03])
        self.assertIsInstance(self.proc1c.Km, tuple)
        self.assertAlmostEqual(self.proc1c.Km[0], 1.204428152)
        self.assertAlmostEqual(self.proc1c.Km[1], 1.204428152 / 4)


class TestRegulatedProcessInVolume(unittest.TestCase):
    def setUp(self):
        self.proc1a = RegulatedProcess({'': 0}, {'B': 1},
                                       k=0.5,
                                       regulating_species='X',
                                       alpha=1,
                                       nH=2,
                                       K50=1e-6,
                                       volume=1e-15)

        self.proc1b = RegulatedProcess.from_string("A -> B",
                                                   k=0.035,
                                                   regulating_species='H',
                                                   alpha=2,
                                                   nH=1,
                                                   K50=(1e-6, 1e-7),
                                                   volume=1e-15)

        self.proc1c = RegulatedProcess.from_string("2A -> B",
                                                   k=3e-4,
                                                   regulating_species='H',
                                                   alpha=0,
                                                   nH=2,
                                                   K50=1e-5,
                                                   volume=1e-15)

    def test_micro_k_vals(self):
        self.assertEqual(self.proc1a.k, 301107038.0)
        self.assertAlmostEqual(self.proc1a.K50, 602.214076)

        self.assertEqual(self.proc1b.k, 0.035)
        self.assertIsInstance(self.proc1b.K50, tuple)
        self.assertAlmostEqual(self.proc1b.K50[0], 602.214076)
        self.assertAlmostEqual(self.proc1b.K50[1], 602.214076 / 10)

        self.assertAlmostEqual(self.proc1c.k, 4.98161720152154e-13)
        self.assertIsInstance(self.proc1c.K50, float)
        self.assertAlmostEqual(self.proc1c.K50, 6022.14076)


class TestRegulatedMichaelisMentenProcessInVolume(unittest.TestCase):
    def setUp(self):
        self.proc1a = RegulatedMichaelisMentenProcess.from_string("F -> G",
                                                                  k=(0.05, 0.02),
                                                                  catalyst='E',
                                                                  Km=2e-8,
                                                                  alpha=2,
                                                                  nH=3,
                                                                  K50=[1e-6, 1e-7],
                                                                  volume=1e-16)

    def test_micro_k_vals(self):
        self.assertTupleEqual(self.proc1a.k, (0.05, 0.02))
        self.assertAlmostEqual(self.proc1a.Km, 1.204428152)
        self.assertIsInstance(self.proc1a.K50, list)
        self.assertAlmostEqual(self.proc1a.K50[0], 60.2214076)
        self.assertAlmostEqual(self.proc1a.K50[1], 6.02214076)


if __name__ == '__main__':
    unittest.main()
