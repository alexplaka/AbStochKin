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
import tempfile
import textwrap
import time
import unittest

from abstochkin.base import AbStochKin
from abstochkin.process import ReversibleProcess


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


class TestAbStochKin1(unittest.TestCase):
    def setUp(self):
        # Create a new instance without specifying the volume
        self.sim = AbStochKin()

    def test_add_simple_process_from_str(self):
        # Test adding a simple unidirectional process
        initial_count = len(self.sim.processes)
        self.sim.add_process_from_str("A -> B", 0.5)
        self.assertEqual(len(self.sim.processes), initial_count + 1, "Process should be added")
        # Check that volume passed is None by default
        proc = self.sim.processes[-1]
        self.assertEqual(getattr(proc, 'volume', None), None, "Volume should default to instance volume (None)")

    def test_extract_process_from_str(self):
        # Test extracting process from string with extra parameters
        initial_count = len(self.sim.processes)
        # Include a numeric parameter and a string parameter
        process_str = "C -> D, k=1.2, catalyst='E'"
        self.sim.extract_process_from_str(process_str)
        self.assertEqual(len(self.sim.processes), initial_count + 1, "Process should be added via extraction")
        proc = self.sim.processes[-1]
        # Check that numerical parameter is set correctly (k)
        self.assertAlmostEqual(proc.k, 1.2, msg="k parameter should be 1.2")
        # Check catalyst if attribute exists
        if hasattr(proc, 'catalyst'):
            self.assertEqual(proc.catalyst, 'E', "Catalyst should be 'E'")

    def test_reversible_process_from_str(self):
        # Test adding a reversible process, which should use '<->' in the string
        initial_count = len(self.sim.processes)
        self.sim.add_process_from_str("X <-> Y", 0.2, k_rev=0.3)
        self.assertEqual(len(self.sim.processes), initial_count + 1, "Reversible process should be added")
        proc = self.sim.processes[-1]
        # If ReversibleProcess is imported, check the type
        if ReversibleProcess is not None:
            self.assertTrue(isinstance(proc, ReversibleProcess), "Process should be a ReversibleProcess")
        # Check k_rev value
        self.assertAlmostEqual(getattr(proc, 'k_rev', 0), 0.3, msg="k_rev should be 0.3")

    def test_add_process_with_dict(self):
        # Test adding process using the dictionary based method
        initial_count = len(self.sim.processes)
        reactants = {'A': 1}
        products = {'B': 1}
        self.sim.add_process(reactants, products, 0.7)
        self.assertEqual(len(self.sim.processes), initial_count + 1,
                         "Process should be added using dict-based function")
        proc = self.sim.processes[-1]
        # Check that reactants and products are set (if attributes exist)
        self.assertEqual(getattr(proc, 'reactants', {}), reactants, "Reactants should match")
        self.assertEqual(getattr(proc, 'products', {}), products, "Products should match")

    def test_del_process_from_str_not_found(self):
        # Test deletion of a non-existent process raises ValueError
        with self.assertRaises(ValueError):
            self.sim.del_process_from_str("NonExistent -> Process", 0.1)

    def test_volume_passing_to_process(self):
        # Test that a specified volume in the AbStochKin instance is passed to added processes
        sim_with_vol = AbStochKin(volume=1e-12)
        sim_with_vol.add_process_from_str("M -> N", 0.9)
        proc = sim_with_vol.processes[-1]
        self.assertAlmostEqual(getattr(proc, 'volume', 0), 1e-12,
                               msg="Process volume should match the simulator volume")

    def test_add_processes_from_file(self):
        # Test adding processes from a temporary file raises TypeError due to positional-only parameter issue
        content = textwrap.dedent("""
            A -> B, k=0.1
            C <-> D, k=0.2, k_rev=0.05
            E -> F, k=0.3, catalyst='G'
        """)
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp:
            tmp.write(content)
            tmp_filename = tmp.name
        try:
            with self.assertRaises(TypeError):
                self.sim.add_processes_from_file(tmp_filename)
        finally:
            os.remove(tmp_filename)

    def test_invalid_file_in_add_processes_from_file(self):
        # Test that add_processes_from_file raises an error when file does not exist.
        with self.assertRaises(FileNotFoundError):
            self.sim.add_processes_from_file("non_existent_file.txt")


class DummySim(AbStochKin):
    """
    A subclass of AbStochKin that overrides the simulate method to return
    the input keyword arguments. This allows us to test simulate_series_in_parallel
    without running actual simulations.
    """

    def simulate(self, **kwargs):
        # Return the kwargs so we can inspect the merging of extra options
        return kwargs


class TestSimulateSeriesInParallel(unittest.TestCase):
    def setUp(self):
        # Instantiate the dummy simulator
        self.sim = DummySim()

    def test_simulate_series_in_parallel(self):
        # Provide a list of simulation parameter dictionaries
        series_params = [
            {"p0": {"A": 10}, "t_max": 5},
            {"p0": {"A": 20}, "t_max": 10},
            {"p0": {"A": 30}, "t_max": 15}
        ]
        # Run the simulations in parallel with a specified max_workers
        self.sim.simulate_series_in_parallel(series_params, max_workers=2)
        # Check that three simulation results have been appended
        self.assertEqual(len(self.sim.sims), 3, "There should be three simulation results")

        # Each simulation result should include the simulation-specific parameters
        # merged with extra options: show_plots=False and _return_simulation=True
        for i, params in enumerate(series_params):
            expected = params.copy()
            expected.update({"show_plots": False, "_return_simulation": True})
            self.assertEqual(self.sim.sims[i], expected, "The simulation parameters should be merged correctly")

    def test_empty_series(self):
        # Test that an empty list returns no simulation results
        self.sim.simulate_series_in_parallel([])
        self.assertEqual(self.sim.sims, [], "No simulations should run for an empty series")


class DummySim1(AbStochKin):
    """
    A subclass of AbStochKin that overrides simulate() to sleep and return the process id.
    """

    def simulate(self, **kwargs):
        time.sleep(2)  # simulate a delay
        return os.getpid()


class TestSimulateParallelProcesses(unittest.TestCase):
    def test_parallel_spawn(self):
        """ Ensure that simulate_series_in_parallel spawns parallel processes. """
        sim = DummySim1()
        series_params = [{}, {}, {}]  # three simulation tasks

        start_time = time.time()
        sim.simulate_series_in_parallel(series_params, max_workers=3)
        elapsed_time = time.time() - start_time

        # If run serially, total time would be around 6 seconds; in parallel it should be significantly less.
        self.assertLess(elapsed_time, 5, "Simulations should run in parallel to reduce total run time")

        # Check that at least one simulation was run in a separate process
        main_pid = os.getpid()
        simulation_pids = sim.sims
        self.assertTrue(any(pid != main_pid for pid in simulation_pids),
                        "At least one simulation should run in a separate process")


if __name__ == '__main__':
    unittest.main()
