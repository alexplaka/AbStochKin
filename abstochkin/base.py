""" Base class, AbStochKin, for initializing and storing all data for
performing stochastic simulations using the Agent-based Kinetics
method. A simulation project can be initialized and run as shown in
the examples below.

Example
-------
>>> from abstochkin import AbStochKin
>>> sim = AbStochKin()
>>> sim.add_process_from_str('A -> ', 0.2)  # degradation process
>>> sim.simulate(p0={'A': 100}, t_max=20)
>>> # All data for the above simulation is stored in `sim.sims[0]`.
>>>
>>> # Now set up a new simulation without actually running it.
>>> sim.simulate(p0={'A': 10}, t_max=10, n=50, run=False)
>>> # All data for the new simulation is stored in `sim.sims[1]`.
>>> # The simulation can then be manually run using methods
>>> # documented in the class `Simulation`.

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

import os
import re
from ast import literal_eval
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Literal

from .het_calcs import get_het_processes
from .process import Process, MichaelisMentenProcess, ReversibleProcess, \
    RegulatedProcess, RegulatedMichaelisMentenProcess
from .simulation import Simulation
from .logging_config import logger
from .utils import log_exceptions

logger = logger.getChild(os.path.basename(__file__))


class AbStochKin:
    """ Base class for Agent-based Kinetics (AbStochKin) simulator.

    Attributes
    ----------
    volume : float, default : None, optional
        The volume *in liters* of the compartment in which the processes
        are taking place.
    volume_unit : str, default : 'L', optional
        A string of the volume unit. The default value is 'L' for liters.
    time_unit : str, default : 'sec', optional
        A string of the time unit to be used for describing the kinetics
        of the given processes.
    processes : list
        A list of the processes that the AbStochKin object has.
    het_processes : list
        A list of the processes where population heterogeneity in one
        of the parameters is to be modeled. This list is a subset of
        the `processes` attribute.
    sims : list
        A list of all simulations performed for the given set of processes.
        Each member of the list is an object of the `Simulation` class and
        contains all data for that simulation.
    """

    def __init__(self,
                 volume: float = None,
                 volume_unit: str = 'L',
                 time_unit: str = 'sec'):
        self.volume = volume
        self.volume_unit = volume_unit
        self.time_unit = time_unit

        self.processes = list()
        self.het_processes = list()
        self.sims = list()

    def add_processes_from_file(self, filename: str):
        """ Add a batch of processes from a text file. """
        with open(filename) as f:
            lines = f.readlines()

        for line in lines:
            self.extract_process_from_str(line)

    def extract_process_from_str(self, process_str: str):
        """
        Extract a process and all of its specified parameters from a string.

        This functions parses a string specifying all values and parameters
        needed to define a process. It then creates a Process object
        based on the extracted data.
        """
        process_str = process_str.replace(' ', '').replace('"', '\'')

        proc_str = process_str.split(',')[0]

        # Extract string parameters first
        patt_str_params = r"(\w+)=('[\w\s,]+')"
        str_params = re.findall(patt_str_params, process_str)

        process_str_remain = re.sub(r"\w+='[\w\s,]+'", '', process_str)

        patt_num_params = r"(\w+)=([\[\(,.\s\d\)\]]+)"
        num_params = re.findall(patt_num_params, process_str_remain)

        all_params = list()
        for name, val_str in num_params + str_params:
            while val_str.endswith(','):
                val_str = val_str[:-1]
            all_params.append((name, literal_eval(val_str)))

        self.add_process_from_str(proc_str, **dict(all_params))

    def add_process_from_str(self,
                             process_str: str,
                             /,
                             k: float | int | list[float | int, ...] | tuple[float | int, float | int],
                             **kwargs):
        """
        Add a process by specifying a string: 'reactants -> products'.
        Additional arguments determine if a specialized process
        (such as a reversible, regulated, or Michaelis-Menten process)
        is to be defined.
        """
        kwargs.setdefault('volume', self.volume)

        if '<->' in process_str or 'k_rev' in kwargs:  # reversible process
            self.processes.append(ReversibleProcess.from_string(process_str, k, **kwargs))
        elif 'catalyst' in kwargs and 'Km' in kwargs and 'regulating_species' not in kwargs:
            self.processes.append(MichaelisMentenProcess.from_string(process_str, k, **kwargs))
        elif 'regulating_species' in kwargs and 'alpha' in kwargs and 'K50' in kwargs and 'nH' in kwargs:
            if 'catalyst' in kwargs and 'Km' in kwargs:  # Regulated MM process
                self.processes.append(RegulatedMichaelisMentenProcess.from_string(process_str,
                                                                                  k, **kwargs))
            else:  # Regulated process
                self.processes.append(RegulatedProcess.from_string(process_str, k, **kwargs))
        else:  # simple unidirectional process
            self.processes.append(Process.from_string(process_str, k, **kwargs))

    def add_process(self,
                    /,
                    reactants: dict,
                    products: dict,
                    k: float | int | list[float | int, ...] | tuple[float | int, float | int],
                    **kwargs):
        """
        Add a process by using a dictionary for the reactants and products.
        Additional arguments determine if a specialized process
        (such as a reversible, regulated, or Michaelis-Menten process)
        is to be defined.
        """
        kwargs.setdefault('volume', self.volume)

        if 'k_rev' in kwargs:  # reversible process
            self.processes.append(ReversibleProcess(reactants, products, k, **kwargs))
        elif 'catalyst' in kwargs and 'Km' in kwargs and 'regulating_species' not in kwargs:
            self.processes.append(MichaelisMentenProcess(reactants, products, k, **kwargs))
        elif 'regulating_species' in kwargs and 'alpha' in kwargs and 'K50' in kwargs and 'nH' in kwargs:
            if 'catalyst' in kwargs and 'Km' in kwargs:  # Regulated MM process
                self.processes.append(
                    RegulatedMichaelisMentenProcess(reactants, products, k, **kwargs))
            else:  # Regulated process
                self.processes.append(RegulatedProcess(reactants, products, k, **kwargs))
        else:  # simple unidirectional process
            self.processes.append(Process(reactants, products, k, **kwargs))

    def del_process_from_str(self,
                             process_str: str,
                             /,
                             k: float | int | list[float | int, ...] | tuple[float | int],
                             **kwargs):
        """ Delete a process by specifying a string: 'reactants -> products'. """
        kwargs.setdefault('volume', self.volume)

        try:
            if '<->' in process_str or 'k_rev' in kwargs:  # reversible process
                self.processes.remove(ReversibleProcess.from_string(process_str, k, **kwargs))
            elif 'catalyst' in kwargs and 'Km' in kwargs and 'regulating_species' not in kwargs:
                self.processes.remove(MichaelisMentenProcess.from_string(process_str, k, **kwargs))
            elif 'regulating_species' in kwargs and 'alpha' in kwargs and 'K50' in kwargs and 'nH' in kwargs:
                if 'catalyst' in kwargs and 'Km' in kwargs:  # Regulated MM process
                    self.processes.remove(
                        RegulatedMichaelisMentenProcess.from_string(process_str, k, **kwargs))
                else:  # Regulated process
                    self.processes.remove(RegulatedProcess.from_string(process_str, k, **kwargs))
            else:  # simple unidirectional process
                self.processes.remove(Process.from_string(process_str, k, **kwargs))
        except ValueError:
            logger.error(f"ValueError: Process to be removed ({process_str}) was not found.")
            raise
        else:
            logger.info(f"Removed: {process_str}, k = {k}, kwargs = {kwargs}")

    def del_process(self,
                    /,
                    reactants: dict,
                    products: dict,
                    k: float | int | list[float | int, ...] | tuple[float | int, float | int],
                    **kwargs):
        """ Delete a process by using a dictionary for the reactants and products. """
        kwargs.setdefault('volume', self.volume)

        try:
            if 'k_rev' in kwargs:  # reversible process
                self.processes.remove(ReversibleProcess(reactants, products, k, **kwargs))
            elif 'catalyst' in kwargs and 'Km' in kwargs and 'regulating_species' not in kwargs:
                self.processes.remove(MichaelisMentenProcess(reactants, products, k, **kwargs))
            elif 'regulating_species' in kwargs and 'alpha' in kwargs and 'K50' in kwargs and 'nH' in kwargs:
                if 'catalyst' in kwargs and 'Km' in kwargs:  # Regulated MM process
                    self.processes.remove(
                        RegulatedMichaelisMentenProcess(reactants, products, k, **kwargs))
                else:  # Regulated process
                    self.processes.remove(RegulatedProcess(reactants, products, k, **kwargs))
            else:  # simple unidirectional process
                self.processes.remove(Process(reactants, products, k, **kwargs))
        except ValueError:
            logger.error(f"ValueError: Process to be removed "
                         f"({reactants=} -> {products=}) was not found.")
            raise
        else:
            lhs, rhs = Process(reactants, products, k)._reconstruct_string()
            logger.info(f"Removed: " + " -> ".join([lhs, rhs]) + f", k = {k}, kwargs = {kwargs}")

    def simulate(self,
                 /,
                 p0: dict,
                 t_max: float | int,
                 dt: float = 0.01,
                 n: int = 100,
                 *,
                 processes: list[Process, ...] = None,
                 random_seed: int = 19,
                 solve_odes: bool = True,
                 ode_method: str = 'RK45',
                 run: bool = True,
                 show_plots: bool = True,
                 plot_backend: Literal['matplotlib', 'plotly'] = 'matplotlib',
                 multithreading: bool = True,
                 max_agents_by_species: dict = None,
                 max_agents_multiplier: int = 2,
                 _return_simulation: bool = False):
        """
        Start an AbStochKin simulation by creating an instance of the class
        `Simulation`. The resulting object is appended to the list
        in the class attribute `AbStochKin.sims`.

        Parameters
        ----------
        p0 : dict[str: int]
            Dictionary specifying the initial population sizes of all
            species in the given processes.
        t_max : float or int
            Numerical value of the end of simulated time in the units
            specified in the class attribute `AbStochKin.time_unit`.
        dt : float, default: 0.1, optional
            The duration of the time interval that the simulation's
            algorithm considers. The current implementation only
            supports a fixed time step interval whose value is `dt`.
        n : int, default: 100, optional
            The number of repetitions of the simulation to be performed.
        processes : list of Process objects, default: None, optional
            A list of the processes to simulate. If the processes to
            simulate are different from the ones in the base class attribute
            `AbStochKin.processes`, they should be specified here.
            If `None`, the processes in the base class attribute
            `AbStochKin.processes` will be used.
        random_seed : int, default: 19, optional
            A number used to seed the random number generator.
        solve_odes : bool, default: True, optional
            Specify whether to numerically solve the system of
            ODEs defined from the given set of processes.
        ode_method : str, default: RK45, optional
            Available ODE methods: RK45, RK23, DOP853, Radau, BDF, LSODA.
        run : bool, default: True, optional
            Specify whether to run an AbStochKin simulation.
        show_plots : bool, default: True, optional
            Specify whether to graph the results of the AbStochKin simulation.
        plot_backend : str, default: 'matplotlib', optional
            `Matplotlib` and `Plotly` are currently supported.
        multithreading : bool, default: True, optional
            Specify whether to parallelize the simulation
            using multithreading. If `False`, the ensemble
            of simulations is run sequentially.
        max_agents_by_species : None or dict, default: dict
            Specification of the maximum number of agents that each
            species should have when running the simulation.
            If `None`, that a default approach will
            be taken by the class `Simulation` and the number
            for each species will be automatically determined
            (see method `Simulation._setup_runtime_data()` for details).
            The entries in the dictionary should be
            `species name (string): number (int)`.
        max_agents_multiplier : float or int, default: 2
            This parameter is used to calculate the maximum number of
            agents of each species that the simulation engine allocates
            memory for. This be determined by multiplying the maximum value
            of the ODE time trajectory for this species by the
            multiplier value specified here.
        _return_simulation : bool
            Determines if the `self.simulate` method returns a `Simulation`
            object or appends it to the list `self.sims`.
            Returning a `Simulation` object is needed when calling the method
            `simulate_series_in_parallel`.
        """

        if max_agents_by_species is None:
            max_agents_by_species = dict()

        # Override `self.processes` if `processes` are specified.
        processes = self.processes if processes is None else processes

        self.het_processes = get_het_processes(processes)

        logger.debug("Creating object of class `Simulation`.")
        sim = Simulation(p0,
                         t_max,
                         dt,
                         n,
                         processes=processes,
                         random_state=random_seed,
                         do_solve_ODEs=solve_odes,
                         ODE_method=ode_method,
                         do_run=run,
                         show_graphs=show_plots,
                         graph_backend=plot_backend,
                         use_multithreading=multithreading,
                         max_agents=max_agents_by_species,
                         max_agents_multiplier=max_agents_multiplier,
                         time_unit=self.time_unit)

        if _return_simulation:
            with log_exceptions():
                assert run, "Must run individual simulations if a series of " \
                            "simulations is to be run with multiprocessing."

            # Set un-pickleable objects to None for data serialization to work
            sim.algo_sequence = None
            sim.progress_bar = None

            return sim
        else:
            self.sims.append(sim)

    def simulate_series_in_parallel(self,
                                    series_kwargs: list[dict[str, Any], ...],
                                    *,
                                    max_workers: int = None):
        """
        Perform a series of simulations in parallel by initializing
        separate processes. Each process runs a simulation and appends
        a `Simulation` object in the list `self.sims`.

        Parameters
        ----------
        series_kwargs : list of dict
            A list containing dictionaries of the keyword arguments for
            performing each simulation in the series. The number of elements
            in the list is the number of simulations that will be run.
        max_workers : int, default: None
            The maximum number of processes to be used for performing
            the given series of simulations. If None, then as many worker
            processes will be created as the machine has processors.

        Examples
        --------
        - Run a series of simulations by varying the initial population size of $A$.
        >>>  from abstochkin import AbStochKin
        >>>
        >>>  sim = AbStochKin()
        >>>  sim.add_process_from_str("A -> B", 0.3, catalyst='E', Km=10)
        >>>  simulation_params = [{"p0": {'A': i, 'B': 0, 'E': 10}, "t_max": 10} for i in range(40, 51)]
        >>>  sim.simulate_series_in_parallel(simulation_params)

        - Run a series of simulations by varying the process parameters.
        >>> from abstochkin import AbStochKin
        >>> from abstochkin.process import Process, RegulatedProcess
        >>>
        >>> simulation_params = []
        >>>
        >>> for α in range(2, 5):
        >>>     procs = [
        >>>         Process.from_string(" -> R", k=0.02),
        >>>         RegulatedProcess.from_string(" -> R",
        >>>                                      k=1,
        >>>                                      regulating_species='S',
        >>>                                      K50=1,
        >>>                                      nH=1,
        >>>                                      alpha=α),
        >>>         Process.from_string("R -> ", k=0.02)
        >>>     ]
        >>>
        >>>     simulation_params.append(
        >>>        dict(
        >>>            processes=procs,
        >>>            p0={'S': 10, "R": 0},
        >>>            t_max=50,
        >>>            dt=0.01,
        >>>            n=20
        >>>        )
        >>>     )
        >>>
        >>> my_sim = AbStochKin()
        >>> my_sim.simulate_series_in_parallel(simulation_params)
        """
        extra_opts = {"show_plots": False, "_return_simulation": True}
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self.simulate, **(kwargs | extra_opts)) for kwargs in
                       series_kwargs]
            for future in futures:
                self.sims.append(future.result())

        logger.info(f"Simulation runs in parallel ({max_workers=}): completed.")
