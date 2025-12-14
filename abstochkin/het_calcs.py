""" Some functions for calculating metrics of population heterogeneity. """

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
from collections import Counter

from numpy import log, ndarray

from .process import Process, ReversibleProcess, MichaelisMentenProcess, \
    RegulatedProcess, RegulatedMichaelisMentenProcess
from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))

ProcessClasses = Process | ReversibleProcess | MichaelisMentenProcess | RegulatedProcess | RegulatedMichaelisMentenProcess


def get_het_processes(processes: list[ProcessClasses, ...]) -> list[ProcessClasses, ...]:
    """
    Filter the heterogeneous processes from a given list of processes.
    A process is heterogeneous if any of its parameters are defined as such.
    """
    het_procs = list()

    for proc in processes:
        if isinstance(proc, MichaelisMentenProcess):
            if proc.is_heterogeneous or proc.is_heterogeneous_Km:
                het_procs.append(proc)
        elif isinstance(proc, RegulatedMichaelisMentenProcess):
            if proc.is_heterogeneous or proc.is_heterogeneous_Km:
                het_procs.append(proc)
            else:
                if isinstance(proc.regulating_species, list):  # multiple regulators
                    if sum(proc.is_heterogeneous_K50) > 0:
                        het_procs.append(proc)
                elif isinstance(proc.regulating_species, str):  # only one regulator
                    if proc.is_heterogeneous_K50:
                        het_procs.append(proc)
        elif isinstance(proc, RegulatedProcess):
            if proc.is_heterogeneous:
                het_procs.append(proc)
            else:
                if isinstance(proc.regulating_species, list):  # multiple regulators
                    if sum(proc.is_heterogeneous_K50) > 0:
                        het_procs.append(proc)
                elif isinstance(proc.regulating_species, str):  # only one regulator
                    if proc.is_heterogeneous_K50:
                        het_procs.append(proc)
        elif isinstance(proc, ReversibleProcess):
            if proc.is_heterogeneous or proc.is_heterogeneous_rev:
                het_procs.append(proc)
        else:  # isinstance(proc, Process)
            if proc.is_heterogeneous:
                het_procs.append(proc)

    return het_procs


def richness(arr: list | ndarray) -> int:
    """
    Calculate the species richness, or how many subspecies
    a species population comprises.
    """
    return len(Counter(arr))


def idx_het(arr: list | ndarray) -> float:
    """
    Calculate the Index of Heterogeneity ($\\psi$), defined as the probability
    that two randomly chosen agents (without replacement) from a species
    population belong to different subspecies.

    - A homogeneous population returns 0.
    - A heterogeneous population with two distinct subspecies of equal
      fractional abundance (χ=0.5) *approaches* 0.5 as the population
      size increases.
    - A fully heterogeneous population returns 1.
    """
    ck = Counter(arr)  # counter of k values
    n_tot = len(arr)  # total number of entries in sequence `arr`

    s = 0
    for v in ck.values():
        s += v * (v - 1)

    try:
        psi = 1 - s / (n_tot * (n_tot - 1))
    except ZeroDivisionError:
        psi = 0

    return psi


def info_het(arr: list | ndarray) -> float:
    """
    Information-theoretic measure of population heterogeneity.

    - A homogeneous population returns 0.
    - A heterogeneous population with two distinct subspecies of equal
      fractional abundance (χ=0.5) returns ln(2). Note that this is
      true regardless of the population size.
    - For a fully heterogeneous population, the measure increases
      with population size and has no upper limit.
    """
    ck = Counter(arr)  # counter of k values
    n_tot = len(arr)  # total number of agents

    s = 0

    for k, v in ck.items():
        chi = v / n_tot
        s -= chi * log(chi)

    return s
