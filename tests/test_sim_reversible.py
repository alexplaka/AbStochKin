import unittest

from abstochkin.base import AbStochKin
from abstochkin.process import Process, ReversibleProcess


class TestReversibleProcessSimulation1(unittest.TestCase):
    def setUp(self):
        self.sim1a = AbStochKin()
        self.sim1a.add_process_from_str("A <-> B", 0.3, k_rev=0.1)
        self.sim1a.simulate(p0={'A': 40, 'B': 10}, t_max=20, dt=0.01, n=100,
                            max_agents_by_species={'A': 60, 'B': 60})

        self.sim1b = AbStochKin()
        self.sim1b.add_process_from_str("A + B <--> C", 0.05, k_rev=0.1)
        self.sim1b.simulate(p0={'A': 40, 'B': 20, 'C': 0}, t_max=10, dt=0.01, n=100,
                            max_agents_by_species={'A': 40, 'B': 40, 'C': 40})

    def test_simulations_1(self):
        # First make sure reversible process is split in two processes:
        # the forward and reverse reactions.
        self.assertEqual(len(self.sim1a.processes), 1)
        self.assertEqual(len(self.sim1a.sims[0]._algo_processes), 2)
        self.assertEqual(type(self.sim1a.processes[0]), ReversibleProcess)
        self.assertEqual(type(self.sim1a.sims[0]._algo_processes[0]), Process)
        self.assertEqual(self.sim1a.sims[0]._algo_processes[0].k, 0.3)
        self.assertEqual(type(self.sim1a.sims[0]._algo_processes[1]), Process)
        self.assertEqual(self.sim1a.sims[0]._algo_processes[1].k, 0.1)
        self.assertEqual(len(self.sim1a.sims[0].algo_sequence), 2)
        # Check simulated trajectories
        self.assertEqual(len(self.sim1a.sims[0]._het_processes), 0)
        self.assertGreaterEqual(self.sim1a.sims[0].results['A']['R^2'], 0.9975)
        self.assertGreaterEqual(self.sim1a.sims[0].results['B']['R^2'], 0.9975)

        # First make sure reversible process is split in two processes:
        # the forward and reverse reactions.
        self.assertEqual(len(self.sim1b.processes), 1)
        self.assertEqual(len(self.sim1b.sims[0]._algo_processes), 2)
        self.assertEqual(type(self.sim1b.processes[0]), ReversibleProcess)
        self.assertEqual(type(self.sim1b.sims[0]._algo_processes[0]), Process)
        self.assertEqual(self.sim1b.sims[0]._algo_processes[0].k, 0.05)
        self.assertEqual(type(self.sim1b.sims[0]._algo_processes[1]), Process)
        self.assertEqual(self.sim1b.sims[0]._algo_processes[1].k, 0.1)
        self.assertEqual(len(self.sim1b.sims[0].algo_sequence), 2)
        # Check simulated trajectories
        self.assertEqual(len(self.sim1b.sims[0]._het_processes), 0)
        self.assertGreaterEqual(self.sim1b.sims[0].results['A']['R^2'], 0.995)
        self.assertGreaterEqual(self.sim1b.sims[0].results['B']['R^2'], 0.995)
        self.assertGreaterEqual(self.sim1b.sims[0].results['C']['R^2'], 0.995)


if __name__ == '__main__':
    unittest.main()
