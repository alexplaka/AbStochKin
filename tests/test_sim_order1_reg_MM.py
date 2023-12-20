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
