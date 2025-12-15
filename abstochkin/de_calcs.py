"""
Perform deterministic calculations on a set of processes.
Construct the ordinary differential equations (ODEs)
describing the system and obtain a numerical solution.
Try to get the system's fixed points numerically and symbolically.
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

from numpy import array, mean
from scipy.integrate import solve_ivp
from scipy.optimize import root
from sympy import Add, Mul, Pow, Matrix, sympify, lambdify, symbols, solve
# We want all single-letter and Greek-letter variables to be symbols.
# We can then use the clashing-symbols dictionaries that have been defined
# as private variables in `_clash` (which includes both single and
# multi-letter names that are defined in `sympy.abc`).
from sympy.abc import _clash

from .process import update_all_species, MichaelisMentenProcess, RegulatedProcess, \
    RegulatedMichaelisMentenProcess
from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


class DEcalcs:
    """
    Perform deterministic calculations given the processes specified
    in an AbStochKin simulation object.

    Attributes
    ----------
    p0 : dict[str: int]
        Dictionary specifying the initial population sizes of all
        species in the given processes.
    t_min : float or int
        Numerical value of the start of simulated time in the units
        specified in the class attribute `time_unit`.
    t_max : float or int
        Numerical value of the end of simulated time in the units
        specified in the class attribute `time_unit`.
    processes : list of Process objects
        A list of all the processes that define the system of ODEs.
        Each process is a `Process` object.
    ode_method: str
        Method for the ODE solver to use.
    time_unit: str
        A string of the time unit to be used for describing the
        kinetics of the given processes.
    """

    def __init__(self,
                 p0: dict,
                 t_min: int | float,
                 t_max: int | float,
                 processes: list,
                 ode_method: str,
                 time_unit: str):
        self.p0 = p0
        self.t_min = t_min
        self.t_max = t_max
        self.processes = processes
        self.ode_method = ode_method
        self.time_unit = time_unit

        # *** For computing and storing the solution to the system of ODEs. ***
        self._all_species, self._procs_by_reactant, self._procs_by_product = update_all_species(
            tuple(self.processes))
        self.odes = dict()  # store species-specific ODEs
        self.odes_matrix = None  # store species-specific ODEs in matrix form
        self._odes_lambdified = None  # for numerical integration using SciPy
        self.jacobian_matrix = None  # Jacobian matrix

        self.odes_sol = None  # numerical solution of species trajectories

        self.fixed_pts_num_sol = None  # numerically-obtained fixed points
        self.fixed_pts_num_sol_info = None  # numerical solver info
        self.fixed_pts_sym_sol = None  # symbolically-obtained fixed points

        self.species_with_ode = list()

        self.setup_ODEs()  # Set up the ODEs and compute the Jacobian matrix

    def setup_ODEs(self, agent_based=True):
        """
        Set up the system of ODEs to be used for computing the
        deterministic trajectory of all the species in the given processes.
        The equations consist of sympy objects.

        Parameters
        ----------
        agent_based : bool, default: True
            If True, set up the agent-based (or microscopic) form of ODEs.
            For instance, for the process `2X -> Y`, the ODE for species
            `X` would include an `X(X - 1)` term instead
            of `X^2` (the canonical form). If False, the canonical form
            of the ODEs is constructed.

        Notes
        -----
        The rate constant, `k`, for a given process is taken to be the mean
        of `k`, unless `k` was defined to be normally-distributed,
        in which case `k` is a 2-tuple and `k[0]` is the specified mean.

        Implementation note: Building up the ODE expressions by separately iterating
        over the processes for products and reactants. This is to properly
        handle 0th order processes for product species. For example, for the
        0th order process ' --> X' with rate constant k1, the ODE is dX/dt = k1.
        """
        # Set up ODEs for product species
        for spec, procs in self._procs_by_product.items():
            terms = []
            for proc in procs:
                k_val = proc.k[0] if isinstance(proc.k, tuple) else mean(proc.k)

                if agent_based:
                    temp = [Mul(*[sympify(name, locals=_clash) - i for i in
                                  range(order)]) if name != '' else 1
                            for name, order in proc.reactants.items()]
                else:
                    temp = [Pow(sympify(name, locals=_clash), order) if name != '' else 1
                            for name, order in proc.reactants.items()]

                new_term = Mul(proc.products[spec], k_val, *temp)
                terms.append(new_term * self.get_term_multiplier(proc))

            if spec in self.odes:
                self.odes[spec] += Add(*terms)
            else:
                self.odes[spec] = Add(*terms)

        # Set up ODEs for reactant species
        for spec, procs in self._procs_by_reactant.items():
            terms = []
            for proc in procs:
                k_val = proc.k[0] if isinstance(proc.k, tuple) else mean(proc.k)

                if agent_based:
                    temp = [Mul(*[sympify(name, locals=_clash) - i for i in range(order)]) for
                            name, order in proc.reactants.items()]
                else:
                    temp = [Pow(sympify(name, locals=_clash), order) for name, order in
                            proc.reactants.items()]

                new_term = Mul(proc.reactants[spec], -1 * k_val, *temp)
                terms.append(new_term * self.get_term_multiplier(proc))

            if spec in self.odes:
                self.odes[spec] += Add(*terms)
            else:
                self.odes[spec] = Add(*terms)

        self.species_with_ode = [sp for sp in self.odes.keys()]

        # Handle species whose population remains constant (e.g., a catalyst)
        if len(self._all_species) != len(self.species_with_ode):
            missing_species = self._all_species - set(self.species_with_ode)
            for m in missing_species:
                self.odes[m] = sympify(0)  # dm/dt = 0
                self.species_with_ode.append(m)

        # Set up the ODEs in matrix form and compute the Jacobian matrix
        self.odes_matrix = Matrix([self.odes[sp] for sp in self.odes.keys()])
        self.jacobian_matrix = self.odes_matrix.jacobian(
            [sympify(sp, locals=_clash) for sp in self.odes.keys()]
        )

        # self.get_fixed_pts_symbolically()

    def lambdify_odes(self):
        """ For converting ODE expressions: sympy -> scipy/numpy. """
        return lambdify([sympify(sp, locals=_clash) for sp in self.odes.keys()],
                        list(self.odes.values()))

    @staticmethod
    def get_term_multiplier(proc):
        """
        Generate the multiplicative term (or terms) needed for generating the
        correct algebraic expressions for specialized processes
        (such as Michaelis-Menten and regulated processes).
        """
        if isinstance(proc, MichaelisMentenProcess):
            Km_val = proc.Km[0] if isinstance(proc.Km, tuple) else mean(proc.Km)
            c = symbols(proc.catalyst)
            return c / (Km_val + sympify(proc.reacts_[0], locals=_clash))
        elif issubclass(type(proc), RegulatedProcess):
            if isinstance(proc.regulating_species, list):
                term = 1
                for sp, a, Kval, hc in zip(proc.regulating_species, proc.alpha, proc.K50, proc.nH):
                    K50_val = Kval[0] if isinstance(Kval, tuple) else mean(Kval)
                    ratio = (symbols(sp) / K50_val) ** hc
                    term *= (1 + a * ratio) / (1 + ratio)
            else:
                K50_val = proc.K50[0] if isinstance(proc.K50, tuple) else mean(proc.K50)
                ratio = (symbols(proc.regulating_species) / K50_val) ** proc.nH
                term = (1 + proc.alpha * ratio) / (1 + ratio)

            if isinstance(proc, RegulatedMichaelisMentenProcess):
                Km_val = proc.Km[0] if isinstance(proc.Km, tuple) else mean(proc.Km)
                c = symbols(proc.catalyst)
                term *= c / (Km_val + sympify(proc.reacts_[0], locals=_clash))

            return term
        else:
            return 1

    def solve_ODEs(self):
        """
        Solve system of ordinary differential equations (ODEs).

        Notes
        -----
        Using the solver `scipy.integrate.solve_ivp`,
        whose method can be one of the following:

        - RK45 : Explicit Runge-Kutta method of order 5(4).
        - RK23 : Explicit Runge-Kutta method of order 3(2).
        - DOP853 : Explicit Runge-Kutta method of order 8.
        - Radau : Implicit Runge-Kutta method of the Radau IIA family of order 5
        - BDF : Implicit multistep variable-order (1 to 5) method based on a
                 backward differentiation formula for the derivative approximation.
        - LSODA : Adams/BDF method with automatic stiffness detection and switching.

        Documentation
        -------------
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        """
        if self._odes_lambdified is None:
            self._odes_lambdified = self.lambdify_odes()

        # Converting initial population values to a list to ensure that the order of
        # species-specific ODEs and initial values are correctly ordered.
        p0_list = [self.p0[sp] for sp in self.odes.keys()]
        self.odes_sol = solve_ivp(lambda t, S: self._odes_lambdified(*S),
                                  t_span=[self.t_min, self.t_max],
                                  y0=p0_list,
                                  method=self.ode_method,
                                  t_eval=None,  # Specify points where the solution is desired
                                  dense_output=True)  # Compute a continuous solution

        self.get_fixed_pts_numerically()

    def get_fixed_pts_numerically(self, root_method="hybr"):
        """ Try to find fixed points numerically.
        See the types of solvers (i.e., methods) to use at documentation below.

        Documentation
        -------------
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
        """
        if self._odes_lambdified is None:
            self._odes_lambdified = self.lambdify_odes()

        fps = root(lambda S: self._odes_lambdified(*S),
                   array([self.p0[sp] for sp in self.odes.keys()]),
                   method=root_method)
        # TODO: May need to do this for more initial guesses to get more/all fixed points.

        if fps.success:
            logger.info(f"(Some) Fixed points found numerically: {list(zip(self.odes.keys(), fps.x))}")
            self.fixed_pts_num_sol = {sp: fp for sp, fp in zip(self.odes.keys(), fps.x)}

        self.fixed_pts_num_sol_info = fps

    def get_fixed_pts_symbolically(self):
        """ Try to solve for the fixed points symbolically.
        Note that one may have to substitute the initial population sizes
        to obtain numerical values here."""
        try:
            self.fixed_pts_sym_sol = solve(self.odes_matrix,
                                           [sympify(sp, locals=_clash) for sp in self.odes.keys()],
                                           dict=True)
        except Exception as e:
            logger.error(f"Fixed points could not be obtained symbolically.\n"
                         f"Exception raised: {e}")
            raise
