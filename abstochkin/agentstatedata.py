"""
Class for storing the state of all agents of a certain species during an
AbStochKin simulation.
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
from copy import deepcopy
from dataclasses import dataclass, field

import numpy as np

from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


@dataclass
class AgentStateData:
    """
    Class for storing the state of all agents of a certain species during an
    AbStochKin simulation.

    Attributes
    ----------
    p_init : int
        The initial population size of the species whose data is
        represented in an `AgentStateData` object.
    max_agents : int
        The maximum number of agents for the species whose data
        is represented in an `AgentStateData` object.
    reps : int
        The number of times the AbStochKin algorithm will repeat a simulation.
        This will be the length of the `asv` list.
    fill_state: int

    asv_ini : numpy.ndarray
        Agent-State Vector (asv) is a species-specific 2-row vector to monitor
        agent state according to Markov's property. This array is the initial
        asv, i.e., at `t=0`. The array shape is `(2, max_agents)`
    asv : list of numpy.ndarray
        A list of length `reps` with copies of `asv_ini`. Each simulation run
        uses its corresponding entry in `asv` to monitor the state of all
        agents.
    """

    p_init: int  # initial population size
    max_agents: int  # maximum number of agents represented in asv
    reps: int  # number of times simulation is repeated
    fill_state: int

    asv_ini: np.ndarray = field(init=False, default_factory=lambda: np.array([]))
    asv: list[np.ndarray] = field(init=False, default_factory=lambda: list(np.array([])))

    def __post_init__(self):
        # Set up initial (t=0) agent-state vector (asv):
        self.asv_ini = np.concatenate(
            (
                np.ones(
                    shape=(2, self.p_init),
                    dtype=np.int8
                ),
                np.full(
                    shape=(2, self.max_agents - self.p_init),
                    fill_value=self.fill_state,
                    dtype=np.int8
                )
            ),
            axis=1
        )

        # Set up separate copy of the initial `asv` for each repetition of the
        # algorithm to facilitate parallelization of ensemble runs.
        self.asv = [deepcopy(self.asv_ini) for _ in range(self.reps)]

    def apply_markov_property(self, r: int):
        """
        The future state of the system depends only on its current state.
        This method is called at the end of each time step in an AbStochKin
        simulation. Therefore, the new agent-state vector becomes the
        current state.
        """
        self.asv[r][0, :] = self.asv[r][1, :]

    def cleanup_asv(self):
        """ Empty the contents of the array `asv`. """
        self.asv = list(np.array([]))

    def get_vals_o1(self,
                    r: int,
                    stream: np.random.Generator,
                    p_vals: np.ndarray,
                    state: int = 1):
        """
        Get random values in [0,1) at a given time step for agents of a given
        state. Agents of other states have a value of zero.

        Get probability values at a given time step for agents of a given state.
        Agents of other states have a transition probability of zero.

        Notes
        -----
        Note that only elements of the `asv` that have the same state in the
        previous and current time steps are considered. This is to ensure that
        agents that have already transitioned to a different state in the
        current time step are not reconsidered for a possible transition.
        """
        nonzero_elems = np.all(self.asv[r] == state, axis=0)
        final_rand_nums = stream.random(self.max_agents) * nonzero_elems
        final_p_vals = p_vals * nonzero_elems

        return final_rand_nums, final_p_vals

    def get_vals_o2(self,
                    other,
                    r: int,
                    stream: np.random.Generator,
                    p_vals: np.ndarray,
                    state: int = 1):
        """
        Get random values in [0,1) at a given time step for interactions between
        agents of a given state. Agents of other states have a value of zero.

        Get probability values at a given time step for interactions between
        agents of a given state. Interactions of agents in other states
        have a transition probability of zero.

        Notes
        -----
        Note that only elements of the `asv` that have the same state in the
        previous and current time steps are considered. This is to ensure that
        agents that have already transitioned to a different state in the
        current time step are not reconsidered for a possible transition.
        """
        nonzero_rows = np.all(self.asv[r] == state, axis=0).reshape(-1, 1)
        nonzero_cols = np.all(other.asv[r] == state, axis=0).reshape(1, -1)

        rand_nums = stream.random(size=(self.max_agents, other.max_agents))

        final_rand_nums = rand_nums * nonzero_rows * nonzero_cols
        final_p_vals = p_vals * nonzero_rows * nonzero_cols

        return final_rand_nums, final_p_vals

    def __str__(self):
        return f"Agent-State Vector with \n" \
               f"Initial population size: {self.p_init}\n" \
               f"Maximum number of agents: {self.max_agents}\n" \
               f"Repeat simulation {self.reps} times\n" \
               f"Fill state: {self.fill_state}"
