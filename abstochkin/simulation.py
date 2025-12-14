"""
Perform an Agent-based Kinetics simulation.

This module contains the code for the class `Simulation`,
which, along with the `SimulationMethodsMixin` class,
does everything that is needed to run an
Agent-based Kinetics simulation and store its results.

The class `AgentStateData` is used by a `Simulation`
object to store and handle some of the necessary runtime data.
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
import contextlib
from typing import Literal

import numpy as np
from tqdm import tqdm

from ._simulation_methods import SimulationMethodsMixin
from .de_calcs import DEcalcs
from .graphing import Graph
from .het_calcs import get_het_processes
from .process import update_all_species, MichaelisMentenProcess, \
    RegulatedProcess, RegulatedMichaelisMentenProcess
from .utils import rng_streams
from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


class Simulation(SimulationMethodsMixin):
    """
    Run an Agent-based Kinetics simulation.

    Attributes
    ----------
    p0 : dict[str: int]
        Dictionary specifying the initial population sizes of all
        species in the given processes.
    t_max : float or int
        Numerical value of the end of simulated time in the units
        specified in the class attribute `AbStochKin.time_unit`.
    dt : float
        The duration of the time interval that the simulation's
        algorithm considers. The current implementation only
        supports a fixed time step interval whose value is `dt`.
    n : int
        The number of repetitions of the simulation to be
        performed.
    processes : list
        A list of the processes to simulate.
    random_state : float or int
        A number used to seed the random number generator.
    use_multithreading : bool
        Specify whether to parallelize the simulation
        using multithreading. If `False`, the ensemble
        of simulations is run sequentially.
    max_agents : dict
        Specification of the maximum number of agents that each
        species should have when running the simulation. An
        empty dictionary signifies that a default approach will
        be taken and the number for each species will be
        automatically determined (see method `_setup_runtime_data()`
        for details). The entries in the dictionary should be
        `species name (string): number (int)`.
    max_agents_multiplier : float or int
        Set the number of possible agents for each species to be the
        maximum value of the ODE solution for the species
        times `max_agents_multiplier`.
    time_unit : str
        A string of the time unit to be used for describing the
        kinetics of the given processes.
    """

    def __init__(self,
                 /,
                 p0: dict,
                 t_max: float,
                 dt: float,
                 n: int,
                 *,
                 processes: list,
                 random_state: int,
                 do_solve_ODEs: bool,
                 ODE_method: str,
                 do_run: bool,
                 show_graphs: bool,
                 graph_backend: Literal['matplotlib', 'plotly'],
                 use_multithreading: bool,
                 max_agents: dict,
                 max_agents_multiplier: float | int,
                 time_unit: str):
        """
        The parameters below are not class attributes, but are part of a
        `Simulation` object's initialization to trigger specific actions
        to be automatically performed. Note that these actions can also
        be performed manually by calling the appropriate methods once a
        class object has been instantiated.

        Other Parameters
        ----------------
        do_solve_ODEs : bool
            If `True`, attempt to numerically solve the system
            of ODEs defined from the given set of processes.
            If `False`, do not attempt to solve the ODEs and
            do not run the simulation.
        ODE_method : str
            Method to use when attempting to solve the system
            of ODEs (if `do_solve_ODEs` is `True`).
            Available ODE methods: RK45, RK23, DOP853, Radau,
            BDF, LSODA.
        do_run : bool
            Specify whether to run the AbStochKin simulation.
            If `False`, then a `Simulation` object is created but
            the simulation is not run. A user can then manually
            run it by calling the class method `run_simulation()`.
        show_graphs : bool
            Specify whether to show graphs of the results.
        graph_backend : str
            `Matplotlib` and `Plotly` are currently supported.
        """

        self.p0 = p0  # dictionary of initial population sizes
        self.t_min = 0  # start simulation at time 0 (by assumption)
        self.t_max = t_max  # end simulation at time t_max
        self.dt = dt  # fixed time interval
        self.n = n  # repeat the simulation this many times
        self.processes = processes
        self.random_state = random_state
        self.use_multithreading = use_multithreading
        self.max_agents = max_agents
        self.max_agents_multiplier = max_agents_multiplier
        self.time_unit = time_unit
        self.graph_backend: Literal['matplotlib', 'plotly'] = graph_backend

        self.all_species, self._procs_by_reactant, self._procs_by_product = \
            update_all_species(tuple(self.processes))

        """ Generate the list of processes to be used by the algorithm. 
        This is done because some processes (e.g., reversible processes) 
        need to be split up into the forward and reverse processes. """
        self._algo_processes = list()
        self._gen_algo_processes(self.processes)
        self._het_processes = get_het_processes(self._algo_processes)
        self._het_processes_num = len(self._het_processes)

        self._validate_p0()  # validate initial population sizes

        # Generate streams of random numbers given a seed
        self.streams = rng_streams(self.n, self.random_state)

        # ******************** Deterministic calculations ********************
        self.de_calcs = DEcalcs(self.p0, self.t_min, self.t_max, self.processes,
                                ode_method=ODE_method, time_unit=self.time_unit)

        if do_solve_ODEs:
            try:
                self.de_calcs.solve_ODEs()
            except Exception as exc:
                logger.exception(f"ODE solver exception:\n{exc}")
                raise
        else:
            if do_run:
                logger.warning("Must specify the maximum number of agents for "
                               "each species when not solving the system ODEs.")
        # ********************************************************************

        self.graph_backend = graph_backend

        self.total_time = self.t_max - self.t_min
        self.t_steps = int(self.total_time / self.dt)
        self.time = np.linspace(0, self.total_time, self.t_steps + 1)

        # Initialize data structures for results, runtime data (rtd),
        # process-specific k values, transition probabilities
        # for 0th and 1st order processes, metrics of population
        # heterogeneity, sequence of functions th algorithm will
        # execute, and progress bar.
        self.results, self.rtd = dict(), dict()

        self.trans_p = dict()  # Process-specific transition probabilities

        # Process-specific parameters that can exhibit heterogeneity
        self.k_vals = dict()
        self.Km_vals = dict()
        self.K50_vals = dict()

        self.k_het_metrics = dict()
        self.Km_het_metrics = dict()
        self.K50_het_metrics = dict()

        # Runtime-specific attributes
        self.algo_sequence = list()
        self.progress_bar = None

        self.graphs = None  # For storing any graphs that are shown

        if do_run:
            logger.info("Running simulation...")
            self.run_simulation()
            if show_graphs:  # Simulation must first be run before plotting
                self.graphs = self.graph_results()

    def run_simulation(self):
        """
        Run an ensemble of simulations and compute statistics
        of simulation data.
        """
        bar_fmt = f"{{desc}}: {{percentage:3.0f}}% |{{bar}}| " \
                  f"n={{n_fmt}}/{{total_fmt}} " \
                  f"[{{elapsed}}{'' if self.use_multithreading else '<{remaining}'}, " \
                  f"{{rate_fmt}}{{postfix}}]"
        self.progress_bar = tqdm(total=self.n,
                                 ncols=65 if self.use_multithreading else 71,
                                 desc="Progress",
                                 bar_format=bar_fmt,
                                 colour='green')
        # Note on progress bar info: Multi-threading makes calculating
        # the remaining time of the simulation unreliable, so only
        # the elapsed time is shown.

        self._setup_data()  # initialize data
        self._gen_algo_sequence()  # generate sequence of processes for algorithm

        self._parallel_run() if self.use_multithreading else self._sequential_run()
        logger.debug("Simulation runs completed.")

        self._compute_trajectory_stats()  # Get statistics on simulation data

        # self._compute_k_het_stats()  # Get statistics on heterogeneity data (k)
        # self._compute_Km_het_stats()  # Get statistics on heterogeneity data (Km)
        # self._compute_K50_het_stats()  # Get statistics on heterogeneity data (Km)
        self._compute_het_stats()  # Get statistics on heterogeneity data
        logger.debug("Statistics computed.")

        self._post_run_cleanup()  # free up some memory
        self.progress_bar.close()

    def graph_results(self,
                      /,
                      graphs_to_show=None,
                      species_to_show=None):
        """
        Make graphs of the results.

        Parameters
        ----------
        species_to_show : `None` or list of string(s), default: `None`
            If `None`, data for all species are plotted.

        graphs_to_show : `None` or string or list of string(s), default: `None`
            If `None`, all graphs are shown. If a string is given then the
            graph that matches the string is shown. A list of strings shows
            all the graphs specified in the list.
            Graph specifications:

            'avg' : Plot the average trajectories and their
                one-standard deviation envelopes. The ODE
                trajectories are also shown.
            'traj' : Plot the species trajectories of individual
                simulations.
            'ode' : Plot the deterministic species trajectories,
                obtained by numerically solving the ODEs.
            'eta' : Plot the coefficient of variation (CoV) and
                the CoV assuming that all processes a species
                participates in obey Poisson statistics.
            'het' : Plot species- and process-specific metrics
                of heterogeneity (`k` and `ψ`) and their
                one-standard deviation envelopes.
        """
        graphs_to_return = list()

        if graphs_to_show is None:
            graphs_to_show = ['avg', 'het']

        if species_to_show is None:
            species_to_show = self.all_species

        if 'avg' in graphs_to_show:
            graph_avg = Graph(backend=self.graph_backend)
            with contextlib.suppress(AttributeError):
                graph_avg.plot_ODEs(self.de_calcs, species=species_to_show, show_plot=False)
            graph_avg.plot_avg_std(self.time, self.results, species=species_to_show)
            graphs_to_return.append(graph_avg)

        if 'traj' in graphs_to_show:
            graph_traj = Graph(backend=self.graph_backend)
            graph_traj.plot_trajectories(self.time, self.results, species=species_to_show)
            graphs_to_return.append(graph_traj)

        if 'ode' in graphs_to_show:
            graph_ode = Graph(backend=self.graph_backend)
            graph_ode.plot_ODEs(self.de_calcs, species=species_to_show)
            graphs_to_return.append(graph_ode)

        if 'eta' in graphs_to_show:
            graph_eta = Graph(backend=self.graph_backend)
            graph_eta.plot_eta(self.time, self.results, species=species_to_show)
            graphs_to_return.append(graph_eta)

        if 'het' in graphs_to_show:
            graph_het = None

            for proc in self._algo_processes:
                if proc.is_heterogeneous:
                    graph_het = Graph(backend=self.graph_backend)
                    graph_het.plot_het_metrics(self.time,
                                               (str(proc), ''),
                                               self.k_het_metrics[proc])

                if isinstance(proc, (MichaelisMentenProcess, RegulatedMichaelisMentenProcess)):
                    if proc.is_heterogeneous_Km:
                        graph_het = Graph(backend=self.graph_backend)
                        graph_het.plot_het_metrics(self.time,
                                                   (str(proc), f"K_m={proc.Km}"),
                                                   self.Km_het_metrics[proc],
                                                   het_attr='Km')

                if isinstance(proc, (RegulatedProcess, RegulatedMichaelisMentenProcess)):
                    if isinstance(proc.is_heterogeneous_K50, list):  # multiple regulators
                        for i in range(len(proc.regulating_species)):
                            if proc.is_heterogeneous_K50[i]:
                                graph_het = Graph(backend=self.graph_backend)
                                extra_str = f"\\textrm{{{proc.regulation_type[i]} by }}" \
                                            f"{proc.regulating_species[i]}, " \
                                            f"K_{{50}}={proc.K50[i]}" # type: ignore
                                graph_het.plot_het_metrics(self.time,
                                                           (str(proc), extra_str),
                                                           self.K50_het_metrics[proc][i],
                                                           het_attr='K50')
                    else:  # one regulator
                        if proc.is_heterogeneous_K50:
                            graph_het = Graph(backend=self.graph_backend)
                            graph_het.plot_het_metrics(self.time,
                                                       (str(proc), f"K_{{50}}={proc.K50}"),
                                                       self.K50_het_metrics[proc],
                                                       het_attr='K50')

            if graph_het is not None:
                graphs_to_return.append(graph_het)

        return graphs_to_return

    def __repr__(self):
        return f"AbStochKin Simulation(p0={self.p0},\n" \
               f"                 t_min={self.t_min},\n" \
               f"                 t_max={self.t_max},\n" \
               f"                 dt={self.dt},\n" \
               f"                 n={self.n},\n" \
               f"                 random_state={self.random_state})"

    def __str__(self):
        descr = f"AbStochKin Simulation object with\n" \
                f"Processes: {self._algo_processes}\n" \
                f"Number of heterogeneous processes: {len(self._het_processes)}\n" \
                f"Initial population sizes: {self.p0}\n" \
                f"Simulated time interval: {self.t_min} - {self.t_max} {self.time_unit} " \
                f"with fixed time step increment {self.dt} {self.time_unit}\n" \
                f"Simulation runs: {self.n}\n" \
                f"Use multi-threading: {self.use_multithreading}\n" \
                f"Random state seed: {self.random_state}\n"

        if self.de_calcs.odes_sol is not None:
            descr += f"Attempt to solve ODEs with method " \
                     f"{self.de_calcs.ode_method}: " \
                     f"{'Successful' if self.de_calcs.odes_sol.success else 'Failed'}\n"

        return descr
