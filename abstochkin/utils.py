"""
Some utility functions on generating random number streams, measuring a
function's runtime, calculating a statistic measuring the goodness of fit
when comparing time series data, and performing unit conversion of
kinetic parameters.
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
import functools
from contextlib import contextmanager
from time import perf_counter

import numpy as np
from numpy import random
from scipy.constants import N_A

from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


def rng_streams(n: int, random_state: int):
    """
    Generate independent streams of random numbers spawned from the
    same initial seed.

    Parameters
    ----------
    n : int
        number of generators/streams to generate.
    random_state : int, optional
        initial seed from which n new seeds are spawned.

    Returns
    -------
    list[PCG64DXSM Generator objects]
        List of `n` generators of independent streams of random numbers

    Notes
    -----
    See https://numpy.org/doc/stable/reference/random/parallel.html for more info.

    On PCG64DXSM:
        - https://numpy.org/doc/stable/reference/random/upgrading-pcg64.html#upgrading-pcg64
        - https://numpy.org/doc/stable/reference/random/bit_generators/pcg64dxsm.html

    Examples
    --------
    >>> rng = rng_streams(5)  # make 5 random number generators
    >>> a = rng[1].integers(1, 10, 100)
    """
    ss = random.SeedSequence(random_state)
    seeds = ss.spawn(n)

    # return [random.default_rng(seed) for seed in seeds]  # PCG64 generator objects

    # Return PCG64DXSM generator objects
    return [random.Generator(random.PCG64DXSM(seed)) for seed in seeds]


def measure_runtime(fcn):
    """ Decorator for measuring the duration of a function's execution. """
    @functools.wraps(fcn)
    def inner(*args, **kwargs):
        """ This is a wrapper function that measures the execution time of the
        decorated function (`fcn`) and logs it. """
        start = perf_counter()
        fcn(*args, **kwargs)
        stop = perf_counter()
        duration = stop - start
        if duration < 60:
            msg = f"Simulation Runtime: {duration:.3f} sec"
        elif 60 <= duration < 3600:
            msg = f"Simulation Runtime: {duration / 60:.3f} min"
        else:
            msg = f"Simulation Runtime: {duration / 3600:.3f} hr"

        logger.info(msg)

    return inner


def r_squared(actual: np.array, theoretical: np.array) -> float:
    """
    Compute the coefficient of determination, $R^2$.

    In the case of comparing the average AbStochKin-simulated species
    trajectory to its deterministic trajectory. Since the latter is only
    meaningful for a homogeneous population, R² should be
    close to `1` for a simulated homogeneous process.
    For a heterogeneous process, it can be interpreted as how close
    the simulated trajectory is to the deterministic trajectory of a
    *homogeneous* process. In this case, $R^2$ would not be expected
    to be close to $1$ and the importance of looking at this metric
    is questionable.

    Parameters
    ----------
    actual : numpy.array
        Actual data obtained through a simulation.
    theoretical : numpy.array
        Theoretical data to compare the actual data to.

    Returns
    -------
    float
        The coefficient of determination, $R^2$.
    """
    # sst: total sum of squares for simulation avg trajectory
    sst = np.nansum((actual - np.nanmean(actual)) ** 2)
    # ssr: sum of square residuals (data vs deterministic prediction)
    ssr = np.nansum((actual - theoretical) ** 2)

    if sst == 0:
        logger.warning("The total sum of squares is zero. R^2 is undefined.")
        return np.nan
    else:
        return 1 - ssr / sst


def macro_to_micro(macro_val: float | int | list[float | int, ...] | tuple[float | int, float | int],
                   volume: float | int,
                   order: int = 0,
                   *,
                   inverse: bool = False) -> float | list[float, ...] | tuple[float, float]:
    """
    Convert a kinetic parameter value from macroscopic to microscopic form.

    The ABK algorithm uses microscopic kinetic constants, thus necessitating
    the conversion of any molar quantities to their microscopic counterpart.
    For a kinetic parameter, the microscopic form is interpreted as the number
    of transition events per second (or whatever the time unit may be).
    For a molar quantity, its microscopic form is the number of particles in
    the given volume.

    Parameters
    ----------
    macro_val : float or int
        The value of the parameter to be converted, expressed in terms of
        molar quantities.
    volume : float or int
        The volume, in liters, in which the process that the given parameter
        value is a descriptor of.
    order : int, default: 0
        The order of the process whose kinetic parameter is to be converted.
        The default value of 0 is for parameters (such as Km or K50) whose
        units are molarity.
    inverse : bool, default: False
        Perform the inverse of this operation. That is, convert
        from microscopic to macroscopic form.

    Returns
    --------
    float
        A kinetic parameter is returned in units of reciprocal seconds.
        A molar quantity is returned as the number of particles in the
        given volume.

    Notes
    -----
    - A kinetic parameter for a 1st order process will remain unchanged
    because its units are already reciprocal seconds.

    Reference
    ---------
    Plakantonakis, Alex. “Agent-based Kinetics: A Nonspatial Stochastic Method
    for Simulating the Dynamics of Heterogeneous Populations.”
    OSF Preprints, 26 July 2019. Web. Section 2.1.
    (Updated version of the paper in "/docs/Agent-basedKinetics_monograph.pdf")
    """
    try:
        assert volume > 0, "The volume has to be a positive quantity."
        assert order >= 0, "The process order cannot be negative."
    except AssertionError:
        logger.error(f"Assertion(s) failed: {volume=}, {order=}\n"
                     f"The volume has to be a positive quantity. "
                     f"The process order cannot be negative.")
        raise

    denom = (N_A * volume) ** (order - 1) if not inverse else 1 / (N_A * volume) ** (order - 1)

    if isinstance(macro_val, (list, tuple)):
        try:
            assert all([True if val >= 0 else False for val in macro_val]), \
                "The parameter values cannot be negative."
            micro_vals = [val / denom for val in macro_val]
            return micro_vals if isinstance(macro_val, list) else tuple(micro_vals)
        except AssertionError:
            logger.error(f"Assertion failed: The parameter values cannot be negative. "
                         f"{macro_val=}")
            raise
    else:
        try:
            assert macro_val >= 0, "The parameter value cannot be negative."
            return macro_val / denom
        except AssertionError:
            logger.error(f"Assertion failed: The parameter value cannot be negative. "
                         f"{macro_val=}")
            raise


@contextmanager
def log_exceptions():
    """
    A context manager for logging exceptions.

    Use this context manager to log any exceptions that occur within
    the given code block. The exception is re-raised after logging.

    Example
    -------
    with log_exceptions():
        # Code that may raise an exception
    """
    try:
        yield
    except Exception as e:
        logger.error("Exception occurred.", exc_info=e)
        raise
