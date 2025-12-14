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
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from functools import partial

import numpy as np

from .agentstatedata import AgentStateData
from .het_calcs import idx_het
from .process import Process, MichaelisMentenProcess, ReversibleProcess, \
    RegulatedProcess, RegulatedMichaelisMentenProcess
from .utils import r_squared, log_exceptions
from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


class SimulationMethodsMixin:
    """ Mixin class with methods necessary for running an AbStochKin simulation. """

    def _validate_p0(self):
        """ A couple of assertions regarding initial population sizes. """
        with log_exceptions():
            assert len(self.p0) == len(self.all_species), \
                f"Specification of initial population sizes does not match number of species.\n" \
                f"len(p0)={len(self.p0)} , len(all_species)={len(self.all_species)}"
            assert all([True if y0 >= 0 else False for y0 in self.p0.values()]), \
                "An initial population size cannot be negative."

    def _setup_data(self):
        """
        Initialize runtime data and `results` dictionary for all species.
        Also initialize process-specific `k` values (`k_vals`),
        transition probabilities (`trans_p`), and metrics of heterogeneity
        (`k_het_metrics`) dictionaries.
        """
        self._setup_runtime_data()

        for sp in self.all_species:
            self.results[sp] = {
                'N': np.empty((self.t_steps + 1, self.n), dtype=np.uint64),
                'N_avg': np.empty(self.t_steps + 1, dtype=np.float64),
                'N_std': np.empty(self.t_steps + 1, dtype=np.float64),
                'eta': np.empty(self.t_steps + 1, dtype=np.float64),  # CoV
                'eta_p': np.empty(self.t_steps + 1, dtype=np.float64),  # Poisson CoV
                'R^2': None
            }
            self.results[sp]['N'][0] = self.p0[sp]

        """ Construct the time-independent PHM.
            Calculate k and transition probability values (done just once). """
        for proc in self._algo_processes:
            if proc.order == 0:
                self.trans_p[proc] = self._get_o0_trans_p(proc)
            elif proc.order == 1:
                self.k_vals[proc], self.trans_p[proc] = self._init_o1_vals(proc)

                if isinstance(proc, (MichaelisMentenProcess, RegulatedMichaelisMentenProcess)):
                    self.Km_vals[proc] = self._init_o1mm_Km_vals(proc)
                if isinstance(proc, (RegulatedProcess, RegulatedMichaelisMentenProcess)):
                    self.K50_vals[proc] = self._init_o1reg_K50_vals(proc)
            else:  # order == 2
                self.k_vals[proc], self.trans_p[proc] = self._init_o2_vals(proc)

                if isinstance(proc, RegulatedProcess):
                    self.K50_vals[proc] = self._init_o2reg_K50_vals(proc)

        self._init_het_metrics()

    def _setup_runtime_data(self):
        """ Initialize runtime data (`rtd`) dictionaries for all species. """

        """ The following is an arbitrary multiplier to account for 
            population size fluctuations and should be sufficient for most 
            cases. The multiplier may need to be adjusted in some cases. """
        if len(self.max_agents) == 0:
            """ Strategy 1: Set the number of possible agents for each 
                species based on the maximum value of the ODE solution for 
                the species times the multiplier. """
            for i, sp in enumerate(self.de_calcs.odes.keys()):
                if sp in self._procs_by_reactant.keys() or sp in self._procs_by_product.keys():
                    sp_max_agents = np.ceil(
                        self.max_agents_multiplier * np.max(self.de_calcs.odes_sol.y[i]))
                else:  # for species whose population size does not change (eg, a catalyst)
                    sp_max_agents = np.ceil(np.max(self.de_calcs.odes_sol.y[i]))

                self.rtd[sp] = AgentStateData(self.p0[sp], int(sp_max_agents), self.n,
                                              fill_state=self._get_fill_state(sp))
                self.max_agents[sp] = int(sp_max_agents)

            """ Strategy 2: Find the maximum from the ODE trajectories of all 
                species and set the number of possible agents for all species 
                to that value times the multiplier. This is an arguably 
                simplistic approach that may result in unnecessary memory (and cpu) 
                usage, however it should be sufficient for most simple systems. """
            # num_agents = int(np.ceil(np.max(self.de_calcs.odes_sol.y))) * max_agents_multiplier
            # for sp in self.all_species:
            #     self.rtd[sp] = AgentStateData(self.p0[sp], num_agents, self.n,
            #                                   fill_state=self._get_fill_state(sp))
            #     self.max_agents[sp] = num_agents
        else:
            """ Strategy 3: Let the user specify the number of possible agents 
                for each species. After examination of the deterministic 
                trajectories and underlying dynamics, it may sometimes 
                be preferred to simply specify the max agents for each species. """
            with log_exceptions():
                assert len(self.max_agents) == len(self.all_species), \
                    "Must specify the maximum number of agents for all species."

            for sp in self.all_species:
                self.rtd[sp] = AgentStateData(self.p0[sp], self.max_agents[sp], self.n,
                                              fill_state=self._get_fill_state(sp))

    def _get_fill_state(self, species_name):
        """
        When setting up the initial agent-state vector (`asv`) of a
        species, agents that are part of the initial population
        have a state of 1 (they are 'alive'). The remaining
        members of the `asv` array are filled with state 0, or -1
        if the species is the product of any 0th order processes.

        Returns
        -------
        fill_state
            0 or -1

        Notes
        -----
        There was a need in an earlier version of the algorithm
        to distinguish between agents that have never been created
        through a 0th order process (or born), and those that have
        been previously born but converted to another species.
        The convention we adopted then is to assign the former
        a state of -1 and the latter a state of 0. We are keeping
        this convention here, although it's not strictly necessary
        for the latest version of the algorithm.
        """
        f_state = 0
        if species_name in self._procs_by_product.keys():
            if 0 in [proc.order for proc in self._procs_by_product[species_name]]:
                f_state = -1
        return f_state

    def _init_het_metrics(self):
        """ Initialize heterogeneity metrics as a list of tuples
        to be converted to a dictionary for each process. """
        # Initialize heterogeneity metrics (for k values) as a list of tuples ...
        # to be converted to a dictionary for each process of a species.
        k_metrics_init = [('k_avg', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                          ('k_std', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                          ('psi', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                          ('<k_avg>', np.zeros(self.t_steps + 1, dtype=np.float64)),
                          ('<k_std>', np.zeros(self.t_steps + 1, dtype=np.float64)),
                          ('psi_avg', np.zeros(self.t_steps + 1, dtype=np.float64)),
                          ('psi_std', np.zeros(self.t_steps + 1, dtype=np.float64))]

        Km_metrics_init = [('Km_avg', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                           ('Km_std', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                           ('psi', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
                           ('<Km_avg>', np.zeros(self.t_steps + 1, dtype=np.float64)),
                           ('<Km_std>', np.zeros(self.t_steps + 1, dtype=np.float64)),
                           ('psi_avg', np.zeros(self.t_steps + 1, dtype=np.float64)),
                           ('psi_std', np.zeros(self.t_steps + 1, dtype=np.float64))]

        K50_metrics_init = [
            ('K50_avg', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
            ('K50_std', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
            ('psi', np.zeros((self.t_steps + 1, self.n), dtype=np.float64)),
            ('<K50_avg>', np.zeros(self.t_steps + 1, dtype=np.float64)),
            ('<K50_std>', np.zeros(self.t_steps + 1, dtype=np.float64)),
            ('psi_avg', np.zeros(self.t_steps + 1, dtype=np.float64)),
            ('psi_std', np.zeros(self.t_steps + 1, dtype=np.float64))]

        for procs in self._procs_by_reactant.values():
            # Compute heterogeneity metrics at t=0
            for proc in procs:
                if proc.is_heterogeneous:
                    k_metrics = deepcopy(k_metrics_init)

                    if proc.order == 1:
                        temp_k_vals = self.k_vals[proc] * (
                                self.rtd[proc.reacts_[0]].asv_ini[0, :] > 0)
                    else:  # order == 2
                        nonzero_rows = (self.rtd[proc.reacts_[0]].asv_ini[1, :] > 0).reshape(-1, 1)
                        nonzero_cols = self.rtd[proc.reacts_[1]].asv_ini[1, :] > 0
                        temp_k_vals = self.k_vals[proc] * nonzero_rows * nonzero_cols

                    nonzero_k_vals = temp_k_vals[np.nonzero(temp_k_vals)]

                    if nonzero_k_vals.size != 0:
                        for i, stat_fcn in enumerate([np.mean, np.std, idx_het]):
                            k_metrics[i][1][0, :] = stat_fcn(nonzero_k_vals)

                    self.k_het_metrics[proc] = dict(k_metrics)

                if isinstance(proc, (MichaelisMentenProcess, RegulatedMichaelisMentenProcess)):
                    if proc.is_heterogeneous_Km:
                        Km_metrics = deepcopy(Km_metrics_init)

                        # Only 1st order Michaelis-Menten processes are currently defined
                        temp_Km_vals = self.Km_vals[proc] * (
                                self.rtd[proc.reacts_[0]].asv_ini[0, :] > 0)

                        nonzero_Km_vals = temp_Km_vals[np.nonzero(temp_Km_vals)]

                        if nonzero_Km_vals.size != 0:
                            for i, stat_fcn in enumerate([np.mean, np.std, idx_het]):
                                Km_metrics[i][1][0, :] = stat_fcn(nonzero_Km_vals)

                        self.Km_het_metrics[proc] = dict(Km_metrics)

                if isinstance(proc, (RegulatedProcess, RegulatedMichaelisMentenProcess)):
                    if isinstance(proc.regulating_species, list):  # multiple regulating species
                        self.K50_het_metrics[proc] = list()

                        for w, Kval in enumerate(proc.K50):
                            if proc.is_heterogeneous_K50[w]:
                                K50_metrics = deepcopy(K50_metrics_init)
                                if proc.order == 1:
                                    temp_K50_vals = self.K50_vals[proc][w] * (
                                            self.rtd[proc.reacts_[0]].asv_ini[0, :] > 0)
                                else:  # order == 2
                                    nonzero_rows = (self.rtd[proc.reacts_[0]].asv_ini[1, :] > 0).reshape(-1, 1)
                                    nonzero_cols = self.rtd[proc.reacts_[1]].asv_ini[1, :] > 0
                                    temp_K50_vals = self.K50_vals[proc][w] * nonzero_rows * nonzero_cols

                                nonzero_K50_vals = temp_K50_vals[np.nonzero(temp_K50_vals)]

                                if nonzero_K50_vals.size != 0:
                                    for i, stat_fcn in enumerate([np.mean, np.std, idx_het]):
                                        K50_metrics[i][1][0, :] = stat_fcn(nonzero_K50_vals)

                                self.K50_het_metrics[proc].append(dict(K50_metrics))
                            else:
                                self.K50_het_metrics[proc].append(None)

                    else:  # only one regulating species
                        if proc.is_heterogeneous_K50:
                            K50_metrics = deepcopy(K50_metrics_init)

                            if proc.order == 1:
                                temp_K50_vals = self.K50_vals[proc] * (
                                        self.rtd[proc.reacts_[0]].asv_ini[0, :] > 0)
                            else:  # order == 2
                                nonzero_rows = (self.rtd[proc.reacts_[0]].asv_ini[1, :] > 0).reshape(-1, 1)
                                nonzero_cols = self.rtd[proc.reacts_[1]].asv_ini[1, :] > 0
                                temp_K50_vals = self.K50_vals[proc] * nonzero_rows * nonzero_cols

                            nonzero_K50_vals = temp_K50_vals[np.nonzero(temp_K50_vals)]

                            if nonzero_K50_vals.size != 0:
                                for i, stat_fcn in enumerate([np.mean, np.std, idx_het]):
                                    K50_metrics[i][1][0, :] = stat_fcn(nonzero_K50_vals)

                            self.K50_het_metrics[proc] = dict(K50_metrics)

    def _get_o0_trans_p(self, proc: Process | RegulatedProcess):
        """ Get the transition probability for a 0th order process. """
        return proc.k * self.dt

    def _init_o1_vals(self, proc: Process):
        """
        Initialize the arrays for the rate constant and
        transition probability values for a 1st order process.
        The latter array is calculated from the former.
        The contents of the arrays depend on whether the process
        is homogeneous or heterogeneous.

        - Homogeneous: The arrays consist of a single repeated value.
        - Heterogeneous: If the population is stratified with a distinct
          number of subspecies, then the array contents are partitioned
          in equal parts according to that number. Each part has the
          value that corresponds to each subspecies. Note that the
          order of the array contents is randomized.
          On the other hand, if the population is normally distributed
          with respect to its rate constant value, then this
          distribution is reflected in the array contents.

        Notes
        -----
        Note that the **sample** mean and standard deviation of the
        generated values in the array will not be exactly the same as
        specified in `proc.k` in the case of normally-distributed
        heterogeneous populations.

        Parameters
        ----------
        proc: Process
            Object of class `Process` whose attribute `order`
            is equal to `1` (i.e., a 1st order process).

        Returns
        -------
        tuple
            k_vals : numpy array
                Contains the values of the rate constant for
                `num_agents` number of agents in a 1st order process.
            p_vals: numpy array
                Contains the values of the transition probabilities
                given a time interval `dt` for a 1st order process.
        """
        num_agents = self.max_agents[proc.reacts_[0]]
        if not proc.is_heterogeneous:  # homogeneous process: one k value
            k_vals = np.full(num_agents, proc.k)
        else:  # heterogeneous process
            if isinstance(proc.k, list):  # distinct subspecies
                k_vals = self._gen_o1_distinct_subspecies_vals(proc, num_agents)
            else:  # proc.k is a tuple: mean, std of normal distribution
                k_vals = self._gen_o1_normal_vals(proc, num_agents)

        p_vals = 1 - np.exp(-1 * k_vals * self.dt)

        return k_vals, p_vals

    def _init_o1mm_Km_vals(self, proc: MichaelisMentenProcess):
        """
        Initialize the array for the Km values for a
        1st order Michaelis-Menten process.
        The contents of the array depends on whether the process
        is homogeneous or heterogeneous.

        - Homogeneous: The array consists of a single repeated value.
        - Heterogeneous: If the population is stratified with a distinct
          number of subspecies, then the array contents are partitioned
          in equal parts according to that number. Each part has the
          value that corresponds to each subspecies. Note that the
          order of the array contents is randomized.
          On the other hand, if the population is normally distributed
          with respect to its Km value, then this
          distribution is reflected in the array contents.
        """
        num_agents = self.max_agents[proc.reacts_[0]]
        if not proc.is_heterogeneous_Km:  # homogeneous process: one Km value
            Km_vals = np.full(num_agents, proc.Km)
        else:  # heterogeneous process wrt Km
            if isinstance(proc.Km, list):  # distinct subspecies wrt Km
                Km_vals = self._gen_o1_distinct_subspecies_vals(proc, num_agents, het_attr='Km')
            else:
                # Assuming Km can only have integer values (output_type).
                Km_vals = self._gen_o1_normal_vals(proc, num_agents, het_attr='Km',
                                                   output_type=int)

        return Km_vals

    def _init_o1reg_K50_vals(self, proc: RegulatedProcess):
        """
        Initialize the array for the K50 values for a regulated
        1st order process.
        The contents of the array depends on whether the process
        is homogeneous or heterogeneous.

        - Homogeneous: The array consists of a single repeated value.
        - Heterogeneous: If the population is stratified with a distinct
          number of subspecies, then the array contents are partitioned
          in equal parts according to that number. Each part has the
          value that corresponds to each subspecies. Note that the
          order of the array contents is randomized.
          On the other hand, if the population is normally distributed
          with respect to its Km value, then this
          distribution is reflected in the array contents.
        """
        num_agents = self.max_agents[proc.reacts_[0]]

        if isinstance(proc.regulating_species, list):  # multiple regulating species
            K50_vals = list()
            for i, Kval in enumerate(proc.K50):
                if not proc.is_heterogeneous_K50[i]:  # homogeneous process: one K50 value
                    K50_vals.append(np.full(num_agents, Kval))
                else:  # heterogeneous process wrt K50
                    if isinstance(Kval, list):  # distinct subspecies wrt K50
                        K50_vals.append(self._gen_o1_distinct_subspecies_vals(proc, num_agents,
                                                                              het_attr='K50',
                                                                              het_attr_idx=i))
                    else:  # K50 is a 2-tuple: normally-distributed
                        # Assuming K50 can only have integer values (output_type).
                        K50_vals.append(self._gen_o1_normal_vals(proc, num_agents,
                                                                 het_attr='K50',
                                                                 het_attr_idx=i,
                                                                 output_type=int))

        else:  # only one regulating species
            if not proc.is_heterogeneous_K50:  # homogeneous process: one K50 value
                K50_vals = np.full(num_agents, proc.K50)
            else:  # heterogeneous process wrt K50
                if isinstance(proc.K50, list):  # distinct subspecies wrt K50
                    K50_vals = self._gen_o1_distinct_subspecies_vals(proc, num_agents,
                                                                     het_attr='K50')
                else:  # K50 is a 2-tuple: normally-distributed
                    # Assuming K50 can only have integer values (output_type).
                    K50_vals = self._gen_o1_normal_vals(proc, num_agents,
                                                        het_attr='K50',
                                                        output_type=int)

        return K50_vals

    def _gen_o1_distinct_subspecies_vals(self,
                                         proc: Process | MichaelisMentenProcess | RegulatedProcess,
                                         num_agents: int,
                                         *,
                                         het_attr: str = 'k',
                                         het_attr_idx: int = None):
        """
        For a heterogeneous population, if the population is stratified
        with a distinct number of subspecies, then the array contents
        are partitioned in equal parts according to that number.
        Each part has the value that corresponds to each subspecies.
        Note that the order of the array contents is randomized.

        Parameters
        ----------
        proc: Process or MichaelisMentenProcess or RegulatedProcess
            Object of class `Process` (or of a class whose parent
            class is `Parent`) whose attribute `order` is equal to
            `1` (i.e., a 1st order process).
        num_agents : int
            The number of elements in the array of values to be
            constructed.
        het_attr : str, default: 'k'
            Specification of the attribute of a process that exhibits
            heterogeneity. The default is the rate constant `k`.
        het_attr_idx : int, default: None
            For cases where the attribute of a process that exhibits
            heterogeneity is a list, this is the index of the element 
            within the list that is desired.

        Returns
        -------
        vals : numpy array
            Contains the values of `het_attr`. The size of the
            array is `num_agents`.
        """
        if het_attr != 'k':
            with log_exceptions():
                assert hasattr(proc, het_attr), f"{proc} does not have attribute {het_attr}."

        attr = getattr(proc, het_attr)
        if isinstance(attr, list) and het_attr_idx is not None:
            attr = attr[het_attr_idx]

        num_subspecies = len(attr)

        with log_exceptions():
            assert num_subspecies <= num_agents, \
                "The number of subspecies cannot be greater than the population size."

        """ `vals` has shape (num_agents,). Note that `num_agents` 
        is the same as `max_agents` when setting up the data. So, this
        sets up the `het_attr` values for all agents that may become 'alive'
        during a simulation. """
        vals = np.array([attr[h] for h in range(num_subspecies) for _ in
                         range(int(num_agents / num_subspecies))])

        """ When the number of agents is not evenly divisible by 
        the number of subspecies, `k_vals` has a shorter length than 
        the required number of agents. That's because 
        `int(num_agents / num_subspecies)` in the above comprehension 
        returns the floor of the ratio. """
        if len(vals) != num_agents:
            num_missing = num_agents - len(vals)
            vals = np.append(vals, [attr[i] for i in range(num_missing)])

        """ Agents in the initial population of the species are 
        assigned `k_vals` sequentially within the array (starting from 
        the beginning). This would over-represent the first (or more) 
        subspecies at the expense of the subsequent ones. 
        Thus, we must shuffle `k_vals` so that the collection of agents 
        in the initial population *approximate* the desired composition 
        of distinct subspecies. """
        np.random.default_rng(seed=self.random_state).shuffle(vals)

        return vals

    def _gen_o1_normal_vals(self,
                            proc: Process | MichaelisMentenProcess | RegulatedProcess,
                            num_agents: int,
                            *,
                            het_attr: str = 'k',
                            het_attr_idx: int = None,
                            output_type: type[float | int] = float,
                            max_tries: int = 1000):
        """ 
        Generate an array of values sampled from a normal distribution.
        Make sure there are no negative values in generated array. 
        
        Parameters
        ----------
        *
        proc: Process or MichaelisMentenProcess or RegulatedProcess
            Object of class `Process` whose attribute `order`
            is equal to `1` (i.e., a 1st order process).
        num_agents: int
            Defines the length of the desired array.
            Effectively represents the number of agents that
            *could* participate in process `proc` during the
            simulation.
        het_attr : str, default: 'k'
            Specification of the process attribute that is 
            normally-distributed. The default, `k`, is the 
            rate constant of the process.
        het_attr_idx : int, default: None
            For cases where the attribute of a process that exhibits
            heterogeneity is a list, this is the index of the element 
            within the list that is desired.
        output_type : type, default: float
            The type of the numbers that will populate the output array.
            The default is `float`, but `int` may sometimes be desired.
        max_tries : int, default: 1000
            Maximum number of tries to generate the array
            without any negative values. If this number is
            exceeded, an exception is raised.
        
        Returns
        -------
        vals : numpy array
            Contains the normally-distributed values of `het_attr` for
            `num_agents` number of agents in a 1st order process. 
            If the array could not be generated after the maximum number 
            of tries, then an exception is raised.
        
        """
        if het_attr != 'k':
            with log_exceptions():
                assert hasattr(proc, het_attr), f"{proc} does not have attribute {het_attr}."

        attr = getattr(proc, het_attr)
        if isinstance(attr, list) and het_attr_idx is not None:
            attr = attr[het_attr_idx]

        rng = np.random.default_rng(seed=self.random_state)

        tries = 0
        while tries <= max_tries:
            vals = rng.normal(*attr, num_agents).astype(output_type)
            if not np.any(vals <= 0):
                return vals
            else:
                tries += 1
        else:
            err_msg = f"Maximum number of tries ({max_tries}) to get non-negative " \
                      f"{het_attr} values in simulated process {proc.__str__()} " \
                      f"has been exceeded. Please reconsider the " \
                      f"distribution of {het_attr} values for this process."
            logger.error(f"Error: {err_msg}")
            raise AssertionError(err_msg)

    def _init_o2_vals(self, proc: Process):
        """
        Initialize the arrays for the rate constant and
        transition probability values for a 2nd order process.
        The latter array is calculated from the former.
        The contents of the arrays depend on whether the process
        is homogeneous or heterogeneous.

        - Homogeneous: The arrays consist of a single repeated value.
        - Heterogeneous: If the population is stratified with a distinct
          number of subspecies, then the array contents are partitioned
          in equal parts according to that number. Each part has the
          value that corresponds to each subspecies. Note that the
          order of the array contents is randomized.
          On the other hand, if the population is normally distributed
          with respect to its rate constant value, then this
          distribution is reflected in the array contents.

        Notes
        -----
        Note that the **sample** mean and standard deviation of the
        generated values in the array will not be exactly the same as
        specified in `proc.k` in the case of normally-distributed
        heterogeneous populations.

        Parameters
        ----------
        proc: Process
            Object of class `Process` whose attribute `order`
            is equal to `2` (i.e., a 2nd order process).

        Returns
        -------
        tuple
            k_vals : numpy array
                Contains the values of the rate constant for
                `num_agents` number of agents in a 1st order process.
            p_vals: numpy array
                Contains the values of the transition probabilities
                given a time interval `dt` for a 1st order process.
        """

        """ phm_shape is a tuple of integers that 
        defines the shape of the population heterogeneity matrix
        (i.e., the `k_vals` and, by extension, `p_vals` arrays).
        Effectively represents the number of agents that
        *could* participate in process `proc` during the
        simulation. For example, for the process A + B -> C,
        the number of rows and columns correspond to the number
        of agents of species A and B respectively. """
        phm_shape = (self.max_agents[proc.reacts_[0]],
                     self.max_agents[proc.reacts_[1]])

        if not proc.is_heterogeneous:  # homogeneous process: one k value
            k_vals = np.full(phm_shape, proc.k)
        else:  # heterogeneous process
            if isinstance(proc.k, list):  # distinct subinteractions
                k_vals = self._gen_o2_distinct_subinteractions_vals(proc, phm_shape)
            else:  # proc.k is a tuple: mean, std of normal distribution
                k_vals = self._gen_o2_normal_vals(proc, phm_shape)

        # If reactants are the same (eg: A + A -> C, a homologous 2nd order
        # process), then make the interaction of an agent with itself impossible.
        if len(set(proc.reacts_)) == 1:
            np.fill_diagonal(k_vals, 0)  # make diagonal entries zero

        p_vals = 1 - np.exp(-1 * k_vals * self.dt)

        return k_vals, p_vals

    def _init_o2reg_K50_vals(self, proc: RegulatedProcess):
        """
        Initialize the array for the K50 values for a regulated
        2nd order process.
        The contents of the array depends on whether the process
        is homogeneous or heterogeneous.

        - Homogeneous: The array consists of a single repeated value.
        - Heterogeneous: If the population is stratified with a distinct
          number of subspecies, then the array contents are partitioned
          in equal parts according to that number. Each part has the
          value that corresponds to each subspecies. Note that the
          order of the array contents is randomized.
          On the other hand, if the population is normally distributed
          with respect to its K50 value, then this
          distribution is reflected in the array contents.
        """
        phm_shape = (self.max_agents[proc.reacts_[0]],
                     self.max_agents[proc.reacts_[1]])

        if isinstance(proc.regulating_species, list):  # multiple regulating species
            K50_vals = list()
            for i, Kval in enumerate(proc.K50):
                if not proc.is_heterogeneous_K50[i]:  # homogeneous process: one K50 value
                    K50_vals.append(np.full(phm_shape, Kval))
                else:  # heterogeneous process wrt K50
                    if isinstance(Kval, list):  # distinct subspecies wrt K50
                        K50_vals.append(self._gen_o2_distinct_subinteractions_vals(proc,
                                                                                   phm_shape,
                                                                                   het_attr='K50',
                                                                                   het_attr_idx=i))
                    else:  # K50 is a 2-tuple: normally-distributed
                        # Assuming K50 can only have integer values (output_type).
                        K50_vals.append(self._gen_o2_normal_vals(proc,
                                                                 phm_shape,
                                                                 het_attr='K50',
                                                                 het_attr_idx=i,
                                                                 output_type=int))

        else:  # only one regulating species
            if not proc.is_heterogeneous_K50:  # homogeneous process: one K50 value
                K50_vals = np.full(phm_shape, proc.K50)
            else:  # heterogeneous process wrt K50
                if isinstance(proc.K50, list):  # distinct subspecies wrt K50
                    K50_vals = self._gen_o2_distinct_subinteractions_vals(proc, phm_shape,
                                                                          het_attr='K50')
                else:  # K50 is a 2-tuple: normally-distributed
                    # Assuming K50 can only have integer values (output_type).
                    K50_vals = self._gen_o2_normal_vals(proc, phm_shape,
                                                        het_attr='K50', output_type=int)

        return K50_vals

    def _gen_o2_distinct_subinteractions_vals(self, proc: Process | RegulatedProcess,
                                              phm_shape: tuple[int, int],
                                              *,
                                              het_attr: str = 'k',
                                              het_attr_idx: int = None):
        """
        For a 2nd order process with distinct subinteractions between
        the reactant agents, as determined by the parameter `het_attr`,
        this function generates a 2-dimensional array of parameter values
        with the desired subinteraction parameter values.

        Parameters
        ----------
        *
        proc : Process or RegulatedProcess
            The 2nd order process whose subinteractions are considered.
        phm_shape : tuple of int, int
            The shape of the population heterogeneity matrix (phm) or
            array to be generated in this function.
        het_attr : str, default: 'k'
            The heterogeneity attribute whose values are to be generated.
            The default is the rate constant `k`. Other possible attributes
            include `K50` for a regulated process.
        het_attr_idx : int, default: None
            For cases where the attribute of a process that exhibits
            heterogeneity is a list, this is the index of the element
            within the list that is desired.

        Notes
        -----
        Note that the **sample** mean and standard deviation of the
        generated values in the array will not be exactly the same as
        expected based on the discrete number of subinteractions.
        """
        if het_attr != 'k':
            with log_exceptions():
                assert hasattr(proc, het_attr), f"{proc} does not have attribute {het_attr}."

        attr = getattr(proc, het_attr)
        if isinstance(attr, list) and het_attr_idx is not None:
            attr = attr[het_attr_idx]

        num_subinteractions = len(attr)

        k_overflow = False
        if len(set(proc.reacts_)) == 1:  # If reactants are the same (eg: A + A -> C)
            if num_subinteractions > phm_shape[0] * phm_shape[1] - np.min(phm_shape):
                k_overflow = True
        else:
            if num_subinteractions > phm_shape[0] * phm_shape[1]:
                k_overflow = True

        if k_overflow:
            error_msg = "The number of subinteractions cannot be greater than the " \
                        "total number of possible inter-agent interactions."
            logger.error(f"Error: {error_msg}")
            raise AssertionError(error_msg)

        """ `k_vals` has shape (num_agents_1, num_agents_2). 
        Note that `num_agents_*` is the same as `max_agents` when setting 
        up the data for a given species. So, this sets up the `het_attr` values 
        for all interactions that may become available during a 
        simulation when agents become 'alive'. """
        rng = np.random.default_rng(seed=self.random_state)
        vals = rng.choice(attr, size=phm_shape, replace=True)
        return vals

    def _gen_o2_normal_vals(self, proc: Process | RegulatedProcess,
                            phm_shape: tuple[int, int],
                            *,
                            het_attr: str = 'k',
                            het_attr_idx: int = None,
                            output_type: type[float | int] = float,
                            max_tries: int = 1000):
        """
        For a 2nd order process where the subinteractions between
        the reactant agents are normally distributed
        with respect to a kinetic parameter value, this function
        generates the 2-dimensional array of parameter values with
        the desired population mean and standard deviation.

        Parameters
        ----------
        proc : Process or RegulatedProcess
            The 2nd order process whose subinteractions are considered.
        phm_shape : tuple of int, int
            The shape of the population heterogeneity matrix (phm) or
            array to be generated in this function.
        het_attr : str, default: 'k'
            The heterogeneity attribute whose values are to be generated.
            The default is the rate constant `k`. Other possible attributes
            include `K50` for a regulated process.
        het_attr_idx : int, default: None
            For cases where the attribute of a process that exhibits
            heterogeneity is a list, this is the index of the element
            within the list that is desired.
        output_type : type
            The type of the numbers that will populate the generated
            array. The default is `float`. For integer values, the
            specified type should be `int`.
        max_tries : int, default: 1000
            The maximum number of times to try generating the array while
            checking that no negative values are in the array given the
            population mean and standard deviation. If this number is
            exceeded, an exception is raised.

        Notes
        -----
        Note that the **sample** mean and standard deviation of the
        generated values in the array will not be exactly the same as
        specified in the `het_attr` process parameter.
        """
        if het_attr != 'k':
            with log_exceptions():
                assert hasattr(proc, het_attr), f"{proc} does not have attribute {het_attr}."

        attr = getattr(proc, het_attr)
        if isinstance(attr, list) and het_attr_idx is not None:
            attr = attr[het_attr_idx]

        rng = np.random.default_rng(seed=self.random_state)

        tries = 0
        while tries <= max_tries:
            vals = rng.normal(*attr, size=phm_shape).astype(output_type)
            if not np.any(vals.flatten() <= 0):
                return vals
            else:
                tries += 1
        else:
            err_msg = f"Maximum number of tries ({max_tries}) to get non-negative " \
                      f"{het_attr} values in simulated process {proc.__str__()} " \
                      f"has been exceeded. Please reconsider the " \
                      f"distribution of {het_attr} values for this process."
            logger.error(f"Error: {err_msg}")
            raise AssertionError(err_msg)

    def _gen_algo_processes(self, procs):
        """
        Generate the list of processes that the algorithm will use
        when running the simulation.

        Notes
        -----
        Convert a reversible process to two separate processes
        representing the forward and reverse reactions. Otherwise,
        keep all irreversible processes unchanged.
        """
        for proc in procs:
            if isinstance(proc, ReversibleProcess):
                forward_proc = Process(proc.reactants, proc.products, proc.k)
                reverse_proc = Process(proc.products, proc.reactants, proc.k_rev)
                self._algo_processes.extend([forward_proc, reverse_proc])
            else:
                self._algo_processes.append(proc)

    def _gen_algo_sequence(self):
        """
        Generate the sequence of functions to be called for the
        algorithm to perform the simulation.
        """
        for proc in self._algo_processes:
            if proc.order == 0:
                if isinstance(proc, RegulatedProcess):
                    if isinstance(proc.regulating_species, list):
                        self.algo_sequence.append(partial(self._order_0_reg_gt1, proc))
                    else:
                        self.algo_sequence.append(partial(self._order_0_reg, proc))
                else:
                    self.algo_sequence.append(partial(self._order_0_base, proc))
            elif proc.order == 1:
                if isinstance(proc, MichaelisMentenProcess):
                    self.algo_sequence.append(partial(self._order_1_mm, proc))
                elif isinstance(proc, RegulatedMichaelisMentenProcess):
                    if isinstance(proc.regulating_species, list):
                        self.algo_sequence.append(partial(self._order_1_reg_mm_gt1, proc))
                    else:
                        self.algo_sequence.append(partial(self._order_1_reg_mm, proc))
                elif isinstance(proc, RegulatedProcess):
                    if isinstance(proc.regulating_species, list):
                        self.algo_sequence.append(partial(self._order_1_reg_gt1, proc))
                    else:
                        self.algo_sequence.append(partial(self._order_1_reg, proc))
                else:
                    self.algo_sequence.append(partial(self._order_1_base, proc))
            else:  # order 2
                if isinstance(proc, RegulatedProcess):
                    if isinstance(proc.regulating_species, list):
                        self.algo_sequence.append(partial(self._order_2_reg_gt1, proc))
                    else:
                        self.algo_sequence.append(partial(self._order_2_reg, proc))
                else:
                    self.algo_sequence.append(partial(self._order_2_base, proc))

    def _parallel_run(self, num_workers=None):
        """
        Run the simulation `n` times in parallel through multithreading.
        Multithreading was chosen instead of multiprocessing because numpy
        functions typically release the GIL while performing their
        computations. In this case, multi-threading allows for parallelism
        without incurring any of the often significant costs in speed related
        to initializing and sharing memory/data between multiple processes.

        From documentation of `concurrent.futures` module at
        https://docs.python.org/3/library/concurrent.futures.html :
        'Default value of max_workers is `min(32, os.cpu_count() + 4)`:
        This default value preserves at least 5 workers for I/O bound tasks.'
        """
        logger.debug("Starting simulation (with multithreading executor)...")
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            executor.map(self._repeat_sim, range(self.n))

    def _sequential_run(self):
        """
        Run the simulation sequentially `n` times
        (without parallelization), i.e., using a simple `for` loop.
        """
        logger.debug("Starting simulation (without multithreading; "
                     "sequential executions of runs)...")
        for i in range(self.n):
            self._repeat_sim(i)

    def _repeat_sim(self, r: int):
        """ Run a repetition (indexed `r`) of the simulation. """
        for t in range(1, self.t_steps + 1):
            for algo_fcn in self.algo_sequence:
                algo_fcn(r, t)

            for sp in self.all_species:
                # Sum up population sizes at end of time step
                self.results[sp]['N'][t, r] = np.sum(self.rtd[sp].asv[r][1, :] > 0)

                # Replace old with new agent-state-vector
                self.rtd[sp].apply_markov_property(r)

        self.progress_bar.update(1)

    def _order_0_base(self, proc: Process, r: int, t: int, transition_p: float = None):
        """
        Assess if a transition event occurs. Then update the agent-state vector
        `asv` accordingly.

        Parameters
        ----------
        proc : Process
            The process object for which the transition is being assessed.
        r : int
            The repetition index of the simulation.
        t : int
            The time step index of the simulation.
        transition_p : float, optional
            The transition probability to use for this event.
            If `None`, defaults to `self.trans_p[proc]`.
        """
        trans_p = transition_p if transition_p is not None else self.trans_p[proc]

        if self.streams[r].random() < trans_p:
            available_agent = np.argmin(self.rtd[proc.prods_[0]].asv[r][0])
            self.rtd[proc.prods_[0]].asv[r][1, available_agent] = 1

    def _order_0_reg(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition probability for a 0th order process
        regulated by a single species
        in a single time step of an AbStochKin simulation.
        """
        ratio = self.results[proc.regulating_species]['N'][t - 1, r] / proc.K50
        mult = (1 + proc.alpha * ratio ** proc.nH) / (1 + ratio ** proc.nH)

        transition_probability = proc.k * mult * self.dt

        self._order_0_base(proc, r, t, transition_probability)

    def _order_0_reg_gt1(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition probability for a regulated 0th order process
        where the number of regulating species in greater than 1,
        in a single time step of an AbStochKin simulation.
        """
        mult = 1
        for i, sp in enumerate(proc.regulating_species):
            ratio = self.results[sp]['N'][t - 1, r] / proc.K50[i]
            mult *= (1 + proc.alpha[i] * ratio ** proc.nH[i]) / (1 + ratio ** proc.nH[i])

        transition_probability = proc.k * mult * self.dt

        self._order_0_base(proc, r, t, transition_probability)

    def _order_1_base(self, proc: Process, r: int, t: int, transition_p: np.ndarray = None):
        """
        Determine the transition events for a 1st order process `proc`
        in a single time step of an AbStochKin simulation. Then update the
        agent-state vector `asv` accordingly.

        Parameters
        ----------
        proc : Process
            The process object for which the transition is being assessed.
        r : int
            The repetition index of the simulation.
        t : int
            The time step index of the simulation.
        transition_p : np.ndarray, optional
            The transition probability to use for this event.
            If `None`, defaults to `self.trans_p[proc]`.
        """
        trans_p = transition_p if transition_p is not None else self.trans_p[proc]

        # Get random numbers (rn) and transition probabilities (tp)
        rn, tp = self.rtd[proc.reacts_[0]].get_vals_o1(r,
                                                       self.streams[r],
                                                       trans_p)

        transition_events = rn < tp  # Determine all transition events in this time step

        # Mark transitions in the ASV of the reactant species
        self.rtd[proc.reacts_[0]].asv[r][1, :] = np.where(transition_events,
                                                          0,
                                                          self.rtd[proc.reacts_[0]].asv[r][1, :])

        num_events = np.sum(transition_events)
        # Mark transitions in the ASV of the product species
        for prod in proc.prods_:
            available_agents = np.nonzero(np.all(self.rtd[prod].asv[r] == 0, axis=0))[0]
            self.rtd[prod].asv[r][1, available_agents[:num_events]] = 1

        if proc.is_heterogeneous:  # Now compute metrics of heterogeneity (k)
            new_k_vals = self.k_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_k_het_metrics(proc, new_k_vals, t, r)

    def _order_1_mm(self, proc: MichaelisMentenProcess, r: int, t: int):
        """
        Determine the transition events for a 1st order Michaelis-Menten
        process `proc` in a single time step of an AbStochKin simulation.
        """
        mult = self.results[proc.catalyst]['N'][t - 1, r] / (
                self.results[proc.reacts_[0]]['N'][t - 1, r] + self.Km_vals[proc])

        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_1_base(proc, r, t, transition_probability)

        if proc.is_heterogeneous_Km:  # Now compute metrics of heterogeneity (Km)
            new_Km_vals = self.Km_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_Km_het_metrics(proc, new_Km_vals, t, r)

    def _order_1_reg(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition events for a 1st order
        process `proc`, regulated by a single species,
        in a single time step of an AbStochKin simulation.
        """
        ratio = self.results[proc.regulating_species]['N'][t - 1, r] / self.K50_vals[proc]
        mult = (1 + proc.alpha * ratio ** proc.nH) / (1 + ratio ** proc.nH)
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_1_base(proc, r, t, transition_probability)

        if proc.is_heterogeneous_K50:  # Now compute metrics of heterogeneity (K50)
            new_K50_vals = self.K50_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_K50_het_metrics(proc, new_K50_vals, t, r)

    def _order_1_reg_gt1(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition events for a 1st order
        process `proc`, regulated by more than 1 species,
        in a single time step of an AbStochKin simulation.
        """
        mult = 1
        for i, sp in enumerate(proc.regulating_species):
            ratio = self.results[sp]['N'][t - 1, r] / self.K50_vals[proc][i]
            mult *= (1 + proc.alpha[i] * ratio ** proc.nH[i]) / (1 + ratio ** proc.nH[i])
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_1_base(proc, r, t, transition_probability)

        for i in range(len(proc.regulating_species)):
            if proc.is_heterogeneous_K50[i]:  # Now compute metrics of heterogeneity (K50)
                new_K50_vals = self.K50_vals[proc][i] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
                self._compute_K50_het_metrics(proc, new_K50_vals, t, r, idx=i)

    def _order_1_reg_mm(self, proc: RegulatedMichaelisMentenProcess, r: int, t: int):
        """
        Determine the transition events for a 1st order
        process `proc`, regulated by a single species and
        obeying Michaelis-Menten kinetics,
        in a single time step of an AbStochKin simulation.
        """
        ratio = self.results[proc.regulating_species]['N'][t - 1, r] / self.K50_vals[proc]
        mult_reg = (1 + proc.alpha * ratio ** proc.nH) / (1 + ratio ** proc.nH)
        mult_mm = self.results[proc.catalyst]['N'][t - 1, r] / (
                self.results[proc.reacts_[0]]['N'][t - 1, r] + self.Km_vals[proc])
        mult = mult_reg * mult_mm
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_1_base(proc, r, t, transition_probability)

        if proc.is_heterogeneous_K50:  # Now compute metrics of heterogeneity (K50)
            new_K50_vals = self.K50_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_K50_het_metrics(proc, new_K50_vals, t, r)

        if proc.is_heterogeneous_Km:  # Now compute metrics of heterogeneity (Km)
            new_Km_vals = self.Km_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_Km_het_metrics(proc, new_Km_vals, t, r)

    def _order_1_reg_mm_gt1(self, proc: RegulatedMichaelisMentenProcess, r: int, t: int):
        """
        Determine the transition events for a 1st order
        process `proc`, regulated by more than one species and
        obeying Michaelis-Menten kinetics,
        in a single time step of an AbStochKin simulation.
        """
        mult_reg = 1
        for i, sp in enumerate(proc.regulating_species):
            ratio = self.results[sp]['N'][t - 1, r] / self.K50_vals[proc][i]
            mult_reg *= (1 + proc.alpha[i] * ratio ** proc.nH[i]) / (1 + ratio ** proc.nH[i])
        mult_mm = self.results[proc.catalyst]['N'][t - 1, r] / (
                self.results[proc.reacts_[0]]['N'][t - 1, r] + self.Km_vals[proc])
        mult = mult_reg * mult_mm
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_1_base(proc, r, t, transition_probability)

        for i in range(len(proc.regulating_species)):
            if proc.is_heterogeneous_K50[i]:  # Now compute metrics of heterogeneity (K50)
                new_K50_vals = self.K50_vals[proc][i] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
                self._compute_K50_het_metrics(proc, new_K50_vals, t, r, idx=i)

        if proc.is_heterogeneous_Km:  # Now compute metrics of heterogeneity (Km)
            new_Km_vals = self.Km_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0)
            self._compute_Km_het_metrics(proc, new_Km_vals, t, r)

    def _order_2_base(self, proc: Process, r: int, t: int, transition_p: np.ndarray = None):
        """
        Determine the transition events for a 2nd order process `proc`
        in a single time step of an AbStochKin simulation. Then update the
        agent-state vector(s) `asv` accordingly.

        Parameters
        ----------
        proc : Process
            The process object for which the transition is being assessed.
        r : int
            The repetition index of the simulation.
        t : int
            The time step index of the simulation.
        transition_p : np.ndarray, optional
            The transition probability to use for this event.
            If `None`, defaults to `self.trans_p[proc]`.
        """
        trans_p = transition_p if transition_p is not None else self.trans_p[proc]

        # Get random numbers (rn) and transition probabilities (tp)
        rn, tp = self.rtd[proc.reacts_[0]].get_vals_o2(self.rtd[proc.reacts_[1]],
                                                       r,
                                                       self.streams[r],
                                                       trans_p)

        transition_events = rn < tp  # Determine all transition events in this time step

        if np.sum(transition_events) > 1:  # IF multiple transition events
            pairs = self._o2_get_unique_pairs(np.argwhere(transition_events))
        else:  # only one or zero transition events
            pairs = np.argwhere(transition_events)

        # Mark transitions in the ASV of each species (or of the same species
        # in the case of a homologous 2nd order process, 2A -> C).
        for pair in pairs:
            self.rtd[proc.reacts_[0]].asv[r][1, pair[0]] = 0
            self.rtd[proc.reacts_[1]].asv[r][1, pair[1]] = 0

        num_events = len(pairs)
        # Mark transitions in the ASV of the product species
        for prod in proc.prods_:
            available_agents = np.nonzero(np.all(self.rtd[prod].asv[r] == 0, axis=0))[0]
            self.rtd[prod].asv[r][1, available_agents[:num_events]] = 1

        # Now compute metrics of heterogeneity for this time step (if applicable)
        if proc.is_heterogeneous:
            new_k_vals = self.k_vals[proc] * (self.rtd[proc.reacts_[0]].asv[r][1, :] > 0).reshape(
                -1, 1) * (self.rtd[proc.reacts_[1]].asv[r][1, :] > 0)
            self._compute_k_het_metrics(proc, new_k_vals, t, r)

    def _order_2_reg(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition events for a 2nd order
        process `proc`, regulated by a single species,
        in a single time step of an AbStochKin simulation.
        """
        ratio = self.results[proc.regulating_species]['N'][t - 1, r] / self.K50_vals[proc]
        mult = (1 + proc.alpha * ratio ** proc.nH) / (1 + ratio ** proc.nH)
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_2_base(proc, r, t, transition_probability)

        if proc.is_heterogeneous_K50:  # Now compute metrics of heterogeneity (K50)
            new_K50_vals = self.K50_vals[proc] * (
                    self.rtd[proc.reacts_[0]].asv[r][1, :] > 0).reshape(-1, 1) * (
                                   self.rtd[proc.reacts_[1]].asv[r][1, :] > 0)
            self._compute_K50_het_metrics(proc, new_K50_vals, t, r)

    def _order_2_reg_gt1(self, proc: RegulatedProcess, r: int, t: int):
        """
        Determine the transition events for a 2nd order
        process `proc`, regulated by more than 1 species,
        in a single time step of an AbStochKin simulation.
        """
        mult = 1
        for i, sp in enumerate(proc.regulating_species):
            ratio = self.results[sp]['N'][t - 1, r] / self.K50_vals[proc][i]
            mult *= (1 + proc.alpha[i] * ratio ** proc.nH[i]) / (1 + ratio ** proc.nH[i])
        transition_probability = 1 - np.exp(-1 * self.k_vals[proc] * mult * self.dt)

        self._order_2_base(proc, r, t, transition_probability)

        for i in range(len(proc.regulating_species)):
            if proc.is_heterogeneous_K50[i]:  # Now compute metrics of heterogeneity (K50)
                new_K50_vals = self.K50_vals[proc][i] * (
                        self.rtd[proc.reacts_[0]].asv[r][1, :] > 0).reshape(-1, 1) * (
                                       self.rtd[proc.reacts_[1]].asv[r][1, :] > 0)
                self._compute_K50_het_metrics(proc, new_K50_vals, t, r, idx=i)

    def _compute_k_het_metrics(self,
                               proc: Process | MichaelisMentenProcess |
                                     RegulatedProcess | RegulatedMichaelisMentenProcess,
                               all_k_vals: np.array,
                               t: int,
                               r: int):
        """
        Calculate metrics of heterogeneity for a 1st or 2nd order process `proc`
        with an array of agent-specific `k` values `all_k_vals`,
        corresponding to time step `t` of repetition `r` of the simulation.
        """
        nonzero_k = all_k_vals[np.nonzero(all_k_vals)]

        if nonzero_k.size != 0:
            self.k_het_metrics[proc]['k_avg'][t, r] = np.mean(nonzero_k)
            self.k_het_metrics[proc]['k_std'][t, r] = np.std(nonzero_k)
            self.k_het_metrics[proc]['psi'][t, r] = idx_het(nonzero_k)

    def _compute_Km_het_metrics(self,
                                proc: MichaelisMentenProcess | RegulatedMichaelisMentenProcess,
                                all_Km_vals: np.array,
                                t: int,
                                r: int):
        """
        Calculate metrics of heterogeneity for a 1st or 2nd order process `proc`
        obeying Michaelis-Menten kinetics with an array of
        agent-specific `Km` values `all_Km_vals`,
        corresponding to time step `t` of repetition `r` of the simulation.
        """
        nonzero_Km = all_Km_vals[np.nonzero(all_Km_vals)]

        if nonzero_Km.size != 0:
            self.Km_het_metrics[proc]['Km_avg'][t, r] = np.mean(nonzero_Km)
            self.Km_het_metrics[proc]['Km_std'][t, r] = np.std(nonzero_Km)
            self.Km_het_metrics[proc]['psi'][t, r] = idx_het(nonzero_Km)

    def _compute_K50_het_metrics(self,
                                 proc: RegulatedProcess | RegulatedMichaelisMentenProcess,
                                 all_K50_vals: np.array,
                                 t: int,
                                 r: int,
                                 *,
                                 idx: int = None):
        """
        Calculate metrics of heterogeneity for a regulated 1st or ***** order process `proc`
        with an array of agent-specific `K50` values `all_K50_vals`,
        corresponding to time step `t` of repetition `r` of the simulation.
        """
        nonzero_K50 = all_K50_vals[np.nonzero(all_K50_vals)]

        if nonzero_K50.size != 0:
            if idx is not None:
                self.K50_het_metrics[proc][idx]['K50_avg'][t, r] = np.mean(nonzero_K50)
                self.K50_het_metrics[proc][idx]['K50_std'][t, r] = np.std(nonzero_K50)
                self.K50_het_metrics[proc][idx]['psi'][t, r] = idx_het(nonzero_K50)
            else:
                self.K50_het_metrics[proc]['K50_avg'][t, r] = np.mean(nonzero_K50)
                self.K50_het_metrics[proc]['K50_std'][t, r] = np.std(nonzero_K50)
                self.K50_het_metrics[proc]['psi'][t, r] = idx_het(nonzero_K50)

    @staticmethod
    def _o2_get_unique_pairs(i_pairs: np.array):
        """
        An agent of a species can only participate in one transition
        event per time step, so this function ensures that multiple
        transition events are not recorded for a given agent.
        For instance, if agent 1 of species A was reported as
        transitioning with agents 5 **and** 9 of species B as partners,
        then only one of those interactions is kept.

        Parameters
        ----------
        i_pairs : numpy.array, shape: (n x 2)
            All interacting pairs of agents that have been reported as
            transitioning in a time step.

        Returns
        -------
        unique_pairs : numpy.array of numpy.array(s) of shape (2, )
            Unique pairs of interacting agents.
        """
        agent_0 = set(i_pairs[:, 0])  # set of number of first agent in pairs
        seen_agent_1 = set()  # set of second agent that is already paired up
        unique_pairs = []

        for a in agent_0:
            for p in i_pairs:
                if p[0] == a:
                    if p[1] not in seen_agent_1:
                        seen_agent_1.add(p[1])  # this agent_1 has been used
                        unique_pairs.append(p)  # add to final list of pairs
                        break
                    else:
                        continue

        return np.array(unique_pairs)

    def _compute_trajectory_stats(self):
        """
        Compute statistics on all the species trajectories
        obtained through an ensemble of `n` simulations.
        """

        for i, sp in enumerate(list(self.de_calcs.odes.keys())):
            self.results[sp]['N_avg'] = np.mean(self.results[sp]['N'], axis=1)
            self.results[sp]['N_std'] = np.std(self.results[sp]['N'], axis=1)

            # Avoid warning about division by zero:
            n_avg_nozeros = np.where(self.results[sp]['N_avg'] == 0,
                                     np.nan, self.results[sp]['N_avg'])

            # Calculate the coefficient of variation, η:
            self.results[sp]['eta'] = self.results[sp]['N_std'] / n_avg_nozeros
            # Calculate η for a Poisson process:
            self.results[sp]['eta_p'] = 1 / np.sqrt(n_avg_nozeros)

            with contextlib.suppress(AttributeError):
                # If the ODEs have been solved, then calculate R^2.
                self.results[sp]['R^2'] = r_squared(self.results[sp]['N_avg'],
                                                    self.de_calcs.odes_sol.sol(self.time).T[:, i])

    def _compute_het_stats(self):
        """ Compute statistics on process-specific metrics of heterogeneity. """
        het_params = [self.k_het_metrics, self.Km_het_metrics, self.K50_het_metrics]
        het_attrs = ['k', 'Km', 'K50']

        for het_param, het_attr in zip(het_params, het_attrs):
            for proc, data in het_param.items():
                if isinstance(data, list):
                    for i, datum in enumerate(data):
                        if datum is not None:
                            het_param[proc][i][f'<{het_attr}_avg>'] = np.mean(datum[f'{het_attr}_avg'], axis=1)
                            het_param[proc][i][f'<{het_attr}_std>'] = np.mean(datum[f'{het_attr}_std'], axis=1)
                            het_param[proc][i]['psi_avg'] = np.mean(datum['psi'], axis=1)
                            het_param[proc][i]['psi_std'] = np.std(datum['psi'], axis=1)
                else:
                    het_param[proc][f'<{het_attr}_avg>'] = np.mean(data[f'{het_attr}_avg'], axis=1)
                    het_param[proc][f'<{het_attr}_std>'] = np.mean(data[f'{het_attr}_std'], axis=1)
                    het_param[proc]['psi_avg'] = np.mean(data['psi'], axis=1)
                    het_param[proc]['psi_std'] = np.std(data['psi'], axis=1)

    def _post_run_cleanup(self):
        """ Delete `asv` for all species after the simulation is done
        to free up the memory associated with this attribute. """
        for sp in self.all_species:
            self.rtd[sp].cleanup_asv()
