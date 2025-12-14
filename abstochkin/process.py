""" Define a process of the form Reactants -> Products. """

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
import re
from typing import Self

from numpy import array

from .utils import macro_to_micro, log_exceptions
from .logging_config import logger

logger = logger.getChild(os.path.basename(__file__))


class Process:
    """
    Define a unidirectional process: Reactants -> Products, where the
    Reactants and Products are specified using standard chemical notation.
    That is, stoichiometric coefficients (integers) and species names are
    specified. For example: 2A + B -> C.

    Attributes
    ----------
    reactants : dict
        The reactants of a given process are specified with
        key-value pairs describing each species name and its
        stoichiometric coefficient, respectively.
    products : dict
        The products of a given process are specified with
        key-value pairs describing each species name and its
        stoichiometric coefficient, respectively.
    k : float, int, list of floats, tuple of floats
        The *microscopic* rate constant(s) for the given process. The data
        type of `k` determines the "structure" of the population as follows:
            - A homogeneous population: if `k` is a single value (float or int),
              then the population is assumed to be homogeneous with all agents
              of the reactant species having kinetics defined by this value.
            - A heterogeneous population with a distinct number of subspecies
              (each with a corresponding `k` value): if `k` is a list of floats,
              then the population is assumed to be heterogeneous with a number
              of subspecies equal to the length of the list.
            - A heterogeneous population with normally-distributed `k` values:
              If `k` is a tuple whose length is 2, then the population is
              assumed to be heterogeneous with a normally distributed `k` value.
              The two entries in the tuple represent the mean and standard
              deviation (in that order) of the desired normal distribution.
    volume : float, default : None, optional
        The volume *in liters* of the compartment in which the processes
        are taking place.
    order : int
        The order of the process (or the molecularity of an elementary process).
        It is the sum of the stoichiometric coefficients of the reactants.
    species : set of strings
        A set of all species in a process.
    reacts_ : list of strings
        A list containing all the reactants in a process.
    prods_ : list of strings
        A list containing all the products in a process.

    Methods
    -------
    from_string
        Class method for creating a Process object from a string.
    """

    def __init__(self,
                 /,
                 reactants: dict[str, int],
                 products: dict[str, int],
                 k: float | int | list[float, ...] | tuple[float, float],
                 *,
                 volume: float | None = None,
                 **kwargs):

        self.reactants = reactants
        self.products = products
        self.k = k
        self.volume = volume

        self._validate_nums()  # make sure there are no errors in given numbers

        # For consistency with processes instantiated using the class method `from_string()`,
        # species denoted as 'None' for a 0th order process are renamed to ''.
        if 'None' in self.reactants.keys():
            self.reactants[''] = self.reactants.pop('None')
        if 'None' in self.products.keys():
            self.products[''] = self.products.pop('None')

        self.order = sum(self.reactants.values())

        self.is_heterogeneous = False if isinstance(self.k, (int, float)) else True

        if self.order == 0:
            with log_exceptions():
                assert not self.is_heterogeneous, \
                    "Since a birth process does not depended on the presence of agents, " \
                    "heterogeneity does not make sense in this context. Please define " \
                    "the rate constant k as a number (not a list or tuple)."

        # Convert macroscopic to microscopic rate constant
        if self.volume is not None:
            self.k = macro_to_micro(self.k, self.volume, self.order)

        # Two ways of storing the involved species:
        # 1) A set of all species
        self.species = set((self.reactants | self.products).keys())
        with contextlib.suppress(KeyError):
            self.species.remove('')  # remove empty species name from any 0th order processes

        # 2) Separate lists
        self.reacts_ = list()  # [reactant species]
        self.prods_ = list()  # [product species]
        self._get_reacts_prods_()

        """ Because Process objects are used as keys in dictionaries used 
        in an AbStochKin simulation, it's much faster to generate the object's 
        string representation once, and then access it whenever it's needed 
        (which could be thousands of times during a simulation). """
        self._str = self.__str__().split(';')[0]

        if len(kwargs) > 0:
            self._lsp(kwargs)

    def _get_reacts_prods_(self):
        """ Make lists of the reactant and product species. Repeated
        elements of a list reflect the order or molecularity of the
        species in the given process. For example, for the process
        `2A + B -> C + D, reacts_ = ['A', 'A', 'B'], prods_ = ['C', 'D']`. """
        for r, m in self.reactants.items():
            for i in range(m):
                self.reacts_.append(r)

        for p, m in self.products.items():
            for i in range(m):
                self.prods_.append(p)

        if '' in self.reacts_:  # remove empty reactant species names
            self.reacts_.remove('')  # from 0th order processes
        if '' in self.prods_:  # remove empty product species names
            self.prods_.remove('')  # from degradation processes

    @classmethod
    def from_string(cls,
                    proc_str: str,
                    /,
                    k: float | int | list[float, ...] | tuple[float, float],
                    *,
                    volume: float | None = None,
                    sep: str = '->',
                    **kwargs) -> Self:
        """ Create a process from a string.

        Parameters
        ----------
        proc_str : str
            A string describing the process in standard chemical notation
            (e.g., 'A + B -> C')
        k : float or int or list of floats or 2-tuple of floats
            The rate constant for the given process. If `k` is a float or
            int, then the process is homogeneous. If `k` is a list, then
            the population of the reactants constsists of distinct subspecies
            or subinteractions depending on the order. If `k` is a 2-tuple,
            then the constant is normally-distributed with a mean and standard
            deviation specified in the tuple's elements.
        volume : float, default : None, optional
            The volume *in liters* of the compartment in which the processes
            are taking place.
        sep : str, default: '->'
            Specifies the characters that distinguish the reactants from the
            products. The default is '->'. The code also treats `-->` as a
            default, if it's present in `proc_str`.

        Notes
        -----
        - Species names should not contain spaces, dashes, and
          should start with a non-numeric character.
        - Zeroth order processes should be specified by an empty space or 'None'.

        Examples
        --------
        >>> Process.from_string("2A + B -> X", 0.3)
        >>> Process.from_string(" -> Y", 0.1)  # for a 0th order (birth) process.
        >>> Process.from_string("Protein_X -> None", 0.15)  # for a 1st order degradation process.
        """

        if len(kwargs) > 0:
            cls._lsp(kwargs)

        sep = '-->' if '-->' in proc_str else sep
        if sep not in proc_str:
            logger.error(f"Cannot find separator '{sep}' in process string: '{proc_str}'")
            raise Exception("Cannot distinguish the reactants from the products.\n"
                            "Please use the *sep* keyword: e.g. sep='->'.")

        lhs, rhs = proc_str.strip().split(sep)  # Left- and Right-hand sides of process
        lhs_terms = lhs.split('+')  # Separate the terms on the left-hand side
        rhs_terms = rhs.split('+')  # Separate the terms on the right-hand side

        return cls(reactants=cls._to_dict(lhs_terms),
                   products=cls._to_dict(rhs_terms),
                   k=k,
                   volume=volume)

    @staticmethod
    def _lsp(kwargs: dict):
        """
        The `Process` class accepts additional arguments (`**kwargs`).
        Since the Process class is a base class for other subclasses,
        this is done so that the Liskov Substitution Principle (LSP)
        is not violated.
        (https://en.wikipedia.org/wiki/Liskov_substitution_principle).
        This way, subclasses override the `from_string` method and have
        additional parameters while remaining consistent with this method
        from their base class. Calling the base instance with the additional
        parameters gives a warning that they will have no effect so that
        the user can intervene, if that's desired.
        """
        msg = f"Warning: Additional parameters {','.join([str(i) for i in kwargs.items()])} " \
              f"will have no effect. "
        if 'k_rev' in kwargs.keys():
            msg += "If that's not what you intended, define the process " \
                   "using ReversibleProcess()."
        if 'regulating_species' in kwargs.keys() or 'alpha' in kwargs.keys() or \
                'nH' in kwargs.keys() or 'K50' in kwargs.keys():
            msg += "If that's not what you intended, define the process " \
                   "using RegulatedProcess()."
        if 'catalyst' in kwargs.keys():
            msg += "If that's not what you intended, define the process " \
                   "using MichaelisMentenProcess()."

        logger.warning(msg)

    @staticmethod
    def _to_dict(terms: list) -> dict:
        """ Convert the information for a side (left, right) of a process
        to a dictionary. """
        side_terms = dict()  # for storing the information of a side of a process
        patt = '^[\\-]*[1-9]+'  # regex pattern (accounts for leading erroneous minus sign)

        if len(terms) == 1 and terms[0].strip().lower() in ['', 'none']:
            spec = ''  # Zeroth order process
            side_terms[spec] = 0
        else:
            for term in terms:
                term = term.strip()
                try:
                    match = re.search(patt, term)
                    stoic_coef = term[slice(*match.span())]  # extract stoichiometric coef
                    spec = re.split(patt, term)[-1].strip()  # extract species name
                    if spec == '' and stoic_coef != 0:
                        logger.error("Cannot find species name.")
                        raise NullSpeciesNameError()
                    stoic_coef = int(stoic_coef)
                except AttributeError:  # when there is no specified stoichiometric coefficient
                    # logger.debug(f"Attribute Error: Cannot find stoichiometric coefficient.\n"
                    #              f"{term=}, Assuming stoichiometric coefficient = 1.")
                    spec = re.split(patt, term)[-1]  # extract species name
                    stoic_coef = 1

                if spec not in side_terms.keys():
                    side_terms[spec] = stoic_coef
                else:
                    side_terms[spec] += stoic_coef

        return side_terms

    def _validate_nums(self):
        """ Make sure coefficients, rate constant, and volume values are not negative. """
        # Check coefficients
        for r, val in (self.reactants | self.products).items():
            with log_exceptions():
                assert val >= 0, f"Coefficient cannot be negative: {val} {r}."

        # Check rate constants
        if isinstance(self.k, (list, tuple)):  # heterogeneous population
            with log_exceptions():
                assert all(array(self.k) >= 0), \
                    f"Rate constant values have to be non-negative: k = {self.k}."
        else:  # when k is a float or int, the population is homogeneous
            with log_exceptions():
                assert self.k >= 0, \
                    f"Rate constant values have to be non-negative: k = {self.k}."

        # For normally-distributed k values, specification is a 2-tuple.
        if isinstance(self.k, tuple):  # normal distribution of k values
            with log_exceptions():
                assert len(self.k) == 2, \
                    "Please specify the mean and standard deviation " \
                    "of k in a 2-tuple: (mean, std)."

        # Check volume
        if self.volume is not None:
            with log_exceptions():
                assert self.volume > 0, f"Volume cannot be negative: {self.volume}."

    def __eq__(self, other):
        if isinstance(other, Process):
            is_equal = (self.k == other.k and
                        self.order == other.order and
                        self.reactants == other.reactants and
                        self.products == other.products and
                        self.species == other.species and
                        self.volume == other.volume)
            return is_equal
        elif isinstance(other, str):
            return self._str == other or self._str.replace(' ', '') == other
        else:
            logger.info(f"{type(self)} and {type(other)} are "
                        f"instances of different classes.")
            return False

    def __hash__(self):
        return hash(self._str)

    def __contains__(self, item):
        return True if item in self.species else False

    def __repr__(self):
        repr_k = macro_to_micro(self.k, self.volume, self.order, inverse=True) if self.volume is not None else self.k
        return f"Process Object: Process.from_string('{self._str.split(',')[0]}', " \
               f"k={repr_k}, " \
               f"volume={self.volume})"

    def __str__(self):
        if isinstance(self.k, (float, int)):
            het_str = "Homogeneous process."
        elif isinstance(self.k, list):
            het_str = f"Heterogeneous process with {len(self.k)} distinct subspecies."
        else:
            het_str = f"Heterogeneous process with normally-distributed k with " \
                      f"mean {self.k[0]} and standard deviation {self.k[1]}."

        lhs, rhs = self._reconstruct_string()

        vol_str = f", volume = {self.volume} L" if self.volume is not None else ""
        return ' -> '.join([lhs, rhs]) + f', k = {self.k}{vol_str}; {het_str}'

    def _reconstruct_string(self):
        lhs = ' + '.join([f"{str(val) + ' ' if val not in [0, 1] else ''}{key}" for key, val in
                          self.reactants.items()])
        rhs = ' + '.join([f"{str(val) + ' ' if val not in [0, 1] else ''}{key}" for key, val in
                          self.products.items()])
        return lhs, rhs


class ReversibleProcess(Process):
    """ Define a reversible process.

    The class-specific attributes are listed below.

    Attributes
    ----------
    k_rev : float or int or list of floats or 2-tuple of floats
        The *microscopic* rate constant for the reverse process.
    is_heterogeneous_rev : bool
        Denotes if the parameter `k_rev` exhibits heterogeneity
        (distinct subspecies/interactions or normally-distributed).

    Notes
    -----
    A `ReversibleProcess` object gets split into two `Process` objects
    (forward and reverse process) when the algorithm runs.
    """

    def __init__(self,
                 /,
                 reactants: dict[str, int],
                 products: dict[str, int],
                 k: float | int | list[float, ...] | tuple[float, float],
                 k_rev: float | int | list[float, ...] | tuple[float, float],
                 *,
                 volume: float | None = None):

        self.k_rev = k_rev  # rate constant for reverse process

        super().__init__(reactants=reactants,
                         products=products,
                         k=k,
                         volume=volume)

        self.is_heterogeneous_rev = False if isinstance(self.k_rev, (int, float)) else True
        self.order_rev = sum(self.products.values())

        if self.volume is not None:  # Convert macroscopic to microscopic rate constants
            self.k_rev = macro_to_micro(k_rev, self.volume, self.order_rev)

    @classmethod
    def from_string(cls,
                    proc_str: str,
                    /,
                    k: float | int | list[float, ...] | tuple[float, float],
                    *,
                    k_rev: float | int | list[float, ...] | tuple[float, float] = 0,
                    volume: float | None = None,
                    sep: str = '<->') -> Self:
        """ Create a reversible process from a string.

        Parameters
        ----------
        proc_str : str
            A string describing the process in standard chemical notation
            (e.g., 'A + B <-> C')
        k : float or int or list of floats or 2-tuple of floats
            The *microscopic* rate constant for the forward process.
        k_rev : float or int or list of floats or 2-tuple of floats
            The *microscopic* rate constant for the reverse process.
        volume : float, default : None, optional
            The volume *in liters* of the compartment in which the processes
            are taking place.
        sep : str, default: '<->'
            Specifies the characters that distinguish the reactants from the
            products. The default is '<->'. The code also treats `<-->` as a
            default, if it's present in `proc_str`.

        Notes
        -----
        - Species names should not contain spaces, dashes, and
          should start with a non-numeric character.

        Examples
        --------
        >>> ReversibleProcess.from_string("2A + B <-> X", 0.3, k_rev=0.2)
        """
        for s in ['<-->', '<=>', '<==>']:
            sep = s if s in proc_str else sep
        if sep not in proc_str:
            logger.error(f"Cannot find separator '{sep}' in process string: '{proc_str}'")
            raise Exception("Cannot distinguish the reactants from the products.\n"
                            "Please use the *sep* keyword: e.g. sep='<->'.")

        lhs, rhs = proc_str.strip().split(sep)  # Left- and Right-hand sides of process
        lhs_terms = lhs.split('+')  # Separate the terms on the left-hand side
        rhs_terms = rhs.split('+')  # Separate the terms on the right-hand side

        return cls(reactants=cls._to_dict(lhs_terms),
                   products=cls._to_dict(rhs_terms),
                   k=k,
                   k_rev=k_rev,
                   volume=volume)

    def __repr__(self):
        repr_k = macro_to_micro(self.k, self.volume, self.order, inverse=True) if self.volume is not None else self.k
        repr_k_rev = macro_to_micro(self.k_rev, self.volume, self.order_rev,
                                    inverse=True) if self.volume is not None else self.k_rev
        return f"ReversibleProcess Object: ReversibleProcess.from_string(" \
               f"'{self._str.split(',')[0]}', " \
               f"k={repr_k}, " \
               f"k_rev={repr_k_rev}, " \
               f"volume={self.volume})"

    def __str__(self):
        if isinstance(self.k, (float, int)):
            het_str = "Forward homogeneous process."
        elif isinstance(self.k, list):
            het_str = f"Forward heterogeneous process with {len(self.k)} " \
                      f"distinct subspecies."
        else:
            het_str = f"Forward heterogeneous process with normally-distributed " \
                      f"k with mean {self.k[0]} and standard deviation {self.k[1]}."

        if isinstance(self.k_rev, (float, int)):
            het_rev_str = "Reverse homogeneous process."
        elif isinstance(self.k_rev, list):
            het_rev_str = f"Reverse heterogeneous process with {len(self.k_rev)} " \
                          f"distinct subspecies."
        else:
            het_rev_str = f"Reverse heterogeneous process with normally-distributed " \
                          f"k with mean {self.k_rev[0]} and standard deviation {self.k_rev[1]}."

        lhs, rhs = self._reconstruct_string()
        vol_str = f", volume = {self.volume} L" if self.volume is not None else ""
        return " <-> ".join([lhs, rhs]) + f", k = {self.k}, k_rev = {self.k_rev}{vol_str}; " \
                                          f"{het_str} {het_rev_str}"

    def _reconstruct_string(self):
        lhs = ' + '.join([f"{str(val) + ' ' if val not in [0, 1] else ''}{key}" for key, val in
                          self.reactants.items()])
        rhs = ' + '.join([f"{str(val) + ' ' if val not in [0, 1] else ''}{key}" for key, val in
                          self.products.items()])
        return lhs, rhs

    def __eq__(self, other):
        if isinstance(other, ReversibleProcess):
            is_equal = (self.k == other.k and
                        self.order == other.order and
                        self.k_rev == other.k_rev and
                        self.order_rev == other.order_rev and
                        self.reactants == other.reactants and
                        self.products == other.products and
                        self.species == other.species and
                        self.volume == other.volume)
            return is_equal
        elif isinstance(other, str):
            return self._str == other or self._str.replace(' ', '') == other
        else:
            logger.info(f"{type(self)} and {type(other)} are instances of different classes.")
            return False

    def __hash__(self):
        return hash(self._str)


class MichaelisMentenProcess(Process):
    """ Define a process that obeys Michaelis-Menten kinetics.

    The class-specific attributes are listed below.

    Attributes
    ----------
    catalyst : str
        Name of the species acting as a catalyst for this process.
    Km : float or int or list of floats or 2-tuple of floats
        *Microscopic* Michaelis constant. Corresponds to the number
        of `catalyst` agents that would produce half-maximal activity.
        Heterogeneity in this parameter is determined by the type of `Km`,
        using the same rules as for parameter `k`.
    is_heterogeneous_Km : bool
        Denotes if the parameter `Km` exhibits heterogeneity
        (distinct subspecies/interactions or normally-distributed).
    """

    def __init__(self,
                 /,
                 reactants: dict[str, int],
                 products: dict[str, int],
                 k: float | int | list[float | int, ...] | tuple[float | int, float | int],
                 *,
                 catalyst: str,
                 Km: float | int | list[float | int, ...] | tuple[float | int, float | int],
                 volume: float | None = None):

        self.catalyst = catalyst
        self.Km = Km

        super().__init__(reactants=reactants,
                         products=products,
                         k=k,
                         volume=volume)

        self.is_heterogeneous_Km = False if isinstance(self.Km, (int, float)) else True
        self.species.add(self.catalyst)
        self._str += f", catalyst = {self.catalyst}, Km = {self.Km}"

        assert self.order != 0, "A 0th order process has no substrate for a catalyst " \
                                "to act on, therefore it cannot follow Michaelis-Menten kinetics."
        if self.order == 2:
            logger.error("2nd order Michaelis-Menten processes are not currently supported.")
            raise NotImplementedError

        if self.volume is not None:  # Convert macroscopic to microscopic Km value
            self.Km = macro_to_micro(Km, self.volume)

    @classmethod
    def from_string(cls,
                    proc_str: str,
                    /,
                    k: float | int | list[float | int, ...] | tuple[float | int, float | int],
                    *,
                    catalyst: str = None,
                    Km: float | int | list[float | int, ...] | tuple[
                        float | int, float | int] = None,
                    volume: float | None = None,
                    sep: str = '->') -> Self:
        """ Create a Michaelis-Menten process from a string.

        Parameters
        ----------
        proc_str : str
            A string describing the process in standard chemical notation
            (e.g., 'A + B -> C')
        k : float or int or list of floats or 2-tuple of floats
            The *microscopic* rate constant for the given process. If `k` is a
            float or int, then the process is homogeneous. If `k` is a list, then
            the population of the reactants constsists of distinct subspecies
            or subinteractions depending on the order. If `k` is a 2-tuple,
            then the constant is normally-distributed with a mean and standard
            deviation specified in the tuple's elements.
        catalyst : str
            Name of species acting as a catalyst.
        Km : float or int or list of floats or 2-tuple of floats
            *Microscopic* Michaelis constant for the process.
            Heterogeneity in this parameter is determined by the type of `Km`,
            using the same rules as for parameter `k`.
        volume : float, default : None, optional
            The volume *in liters* of the compartment in which the processes
            are taking place.
        sep : str, default: '->'
            Specifies the characters that distinguish the reactants from the
            products. The default is '->'. The code also treats `-->` as a
            default, if it's present in `proc_str`.

        Notes
        -----
        - Species names should not contain spaces, dashes, and
          should start with a non-numeric character.
        - Zeroth order processes should be specified by an empty space or 'None'.

        Examples
        --------
        >>> MichaelisMentenProcess.from_string("A -> X", k=0.3, catalyst='E', Km=10)
        >>> MichaelisMentenProcess.from_string("A -> X", k=0.3, catalyst='alpha', Km=(10, 1))
        """
        sep = '-->' if '-->' in proc_str else sep
        if sep not in proc_str:
            logger.error(f"Cannot find separator '{sep}' in process string: '{proc_str}'")
            raise Exception("Cannot distinguish the reactants from the products.\n"
                            "Please use the *sep* keyword: e.g. sep='->'.")

        lhs, rhs = proc_str.strip().split(sep)  # Left- and Right-hand sides of process
        lhs_terms = lhs.split('+')  # Separate the terms on the left-hand side
        rhs_terms = rhs.split('+')  # Separate the terms on the right-hand side

        return cls(reactants=cls._to_dict(lhs_terms),
                   products=cls._to_dict(rhs_terms),
                   k=k,
                   catalyst=catalyst,
                   Km=Km,
                   volume=volume)

    def __repr__(self):
        repr_k = macro_to_micro(self.k, self.volume, self.order, inverse=True) if self.volume is not None else self.k
        repr_Km = macro_to_micro(self.Km, self.volume, inverse=True) if self.volume is not None else self.Km
        return f"MichaelisMentenProcess Object: " \
               f"MichaelisMentenProcess.from_string('{self._str.split(',')[0]}', " \
               f"k={repr_k}, " \
               f"catalyst='{self.catalyst}', " \
               f"Km={repr_Km}, " \
               f"volume={self.volume})"

    def __str__(self):
        if isinstance(self.Km, (float, int)):
            Km_het_str = "Homogeneous process with respect to Km."
        elif isinstance(self.k, list):
            Km_het_str = f"Heterogeneous process with respect to Km " \
                         f"with {len(self.Km)} distinct subspecies."
        else:
            Km_het_str = f"Heterogeneous process with normally-distributed Km with " \
                         f"mean {self.Km[0]} and standard deviation {self.Km[1]}."

        return super().__str__() + f" Catalyst: {self.catalyst}, " \
                                   f"Km = {self.Km}, {Km_het_str}"

    def __eq__(self, other):
        if isinstance(other, MichaelisMentenProcess):
            is_equal = (self.k == other.k and
                        self.order == other.order and
                        self.reactants == other.reactants and
                        self.products == other.products and
                        self.catalyst == other.catalyst and
                        self.Km == other.Km and
                        self.species == other.species and
                        self.volume == other.volume)
            return is_equal
        elif isinstance(other, str):
            return self._str == other or self._str.replace(' ', '') == other
        else:
            logger.info(f"{type(self)} and {type(other)} are instances of different classes.")
            return False

    def __hash__(self):
        return hash(self._str)


class RegulatedProcess(Process):
    """ Define a process that is regulated.

    This class allows a Process to be defined in terms of how it is regulated.
    If there is only one regulating species, then the parameters have the same
    type as would be expected for a homogeneous/heterogeneous process.
    If there are multiple regulating species, then all parameters are a list
    of their expected type, with the length of the list being equal to the
    number of regulating species.

    The class-specific attributes (except for `k`, which requires some
    additional notes) are listed below.
    
    Attributes
    ----------
    k : float or int or list of floats or 2-tuple of floats
        The *microscopic* rate constant for the given process. It is the *basal*
        rate constant in the case of activation (or the minimum `k` value)
        and the maximum rate constant in the case of repression.
    regulating_species : str or list of str
        Name of the regulating species. Multiple species can be specified as
        comma-separated in a string or a list of strings with the species names.
    alpha : float or int or list[float or int]
        Parameter denoting the degree of activation/repression.

            - 0 <= alpha < 1: repression
            - alpha = 1: no regulation
            - alpha > 1: activation
            
        alpha is a multiplier: in the case of activation, the maximum 
        rate constant value will be `alpha * k`. 
        In the case of repression, the minimum 
        rate constant value will be `alpha * k`. 
    K50 : float or int or list of floats or 2-tuple of floats or list[float or int or list of floats or 2-tuple of floats]
        *Microscopic* constant that corresponds to the number of
        `regulating_species` agents that would produce 
        half-maximal activation/repression. 
        Heterogeneity in this parameter is determined by the type of `K50`,
        using the same rules as for parameter `k`.
    nH : float or int or list[float or int]
        Hill coefficient for the given process. Indicates the degree of 
        cooperativity in the regulatory interaction. 
    is_heterogeneous_K50 : bool or list of bool
        Denotes if the parameter `K50` exhibits heterogeneity
        (distinct subspecies/interactions or normally-distributed).
    regulation_type : str or list of str
        The type of regulation for this process based on the value of alpha:
        'activation' or 'repression' or 'no regulation'.

    Notes
    -----
    Allowing a 0th order process to be regulated. However, heterogeneity
    in `k` and `K50` (or any other parameters) is not allowed for such
    a process.
    """

    def __init__(self,
                 /,
                 reactants: dict[str, int],
                 products: dict[str, int],
                 k: float | int | list[float, ...] | tuple[float, float],
                 *,
                 regulating_species: str | list[str, ...],
                 alpha: float | int | list[float | int, ...],
                 K50: float | int | list[float | int, ...] | tuple[float | int, float | int] |
                      list[float | int | list[float | int, ...] | tuple[float | int, float | int]],
                 nH: float | int | list[float | int, ...],
                 volume: float | None = None):

        if isinstance(regulating_species, str):
            reg_sp_list = regulating_species.replace(' ', '').split(',')
            self.regulating_species = reg_sp_list[0] if len(reg_sp_list) == 1 else reg_sp_list
        else:  # if it is a list
            self.regulating_species = regulating_species

        self.alpha = alpha
        self.K50 = K50
        self.nH = nH

        if isinstance(K50, list):
            self.is_heterogeneous_K50 = [False if isinstance(val, (int, float)) else True for val
                                         in K50]
        else:
            self.is_heterogeneous_K50 = False if isinstance(self.K50, (int, float)) else True

        super().__init__(reactants=reactants,
                         products=products,
                         k=k,
                         volume=volume)

        if self.volume is not None:  # Convert macroscopic to microscopic K50 value
            self.K50 = macro_to_micro(K50, self.volume)

        self._str += f", regulating_species = {self.regulating_species}, alpha = {self.alpha}, " \
                     f"K50 = {self.K50}, nH = {self.nH}"

        self._validate_reg_params()

        if isinstance(self.alpha, list):
            self.regulation_type = list()
            for a in self.alpha:
                reg_type = 'activation' if a > 1 else 'repression' if a < 1 else 'no regulation'
                self.regulation_type.append(reg_type)
        else:
            self.regulation_type = 'activation' if self.alpha > 1 else 'repression' if self.alpha < 1 else 'no regulation'

    def _validate_reg_params(self):
        """ Validate the parameters specific to the regulation. """
        if isinstance(self.regulating_species, list):  # multiple regulating species
            # First check that the right number of values for each parameter are specified
            rs_num = len(self.regulating_species)
            msg = f"Must specify {rs_num} values when there are {rs_num} regulating species."
            with log_exceptions():
                assert len(self.alpha) == rs_num, msg.replace('#', 'alpha')
                assert len(self.K50) == rs_num, msg.replace('#', 'K50')
                assert len(self.nH) == rs_num, msg.replace('#', 'nH')

            for i in range(len(self.regulating_species)):
                with log_exceptions():
                    assert self.alpha[i] >= 0, "The alpha parameter must be nonnegative."

                if self.alpha[i] == 1:
                    logger.warning("Warning: alpha=1 means the process is not regulated.")

                if isinstance(self.K50[i], (float, int)):
                    with log_exceptions():
                        assert self.K50[i] > 0, "K50 has to be positive."
                elif isinstance(self.K50[i], list):
                    with log_exceptions():
                        assert all([True if val > 0 else False for val in self.K50[i]]), \
                            "Subspecies K50 values have to be positive."
                else:  # isinstance(self.K50, tuple)
                    with log_exceptions():
                        assert self.K50[i][0] > 0 and self.K50[i][1] > 0, \
                            "Mean and std of K50 have to be positive."

                if self.order == 0:
                    with log_exceptions():
                        assert not self.is_heterogeneous_K50[i], \
                            "Heterogeneity in parameter K50 is not allowed for a 0th order process."

                with log_exceptions():
                    assert self.nH[i] > 0, "nH has to be positive."

        else:  # just one regulating species
            with log_exceptions():
                assert self.alpha >= 0, "The alpha parameter must be nonnegative."

            if self.alpha == 1:
                logger.warning("Warning: alpha=1 means the process is not regulated.")

            if isinstance(self.K50, (float, int)):
                with log_exceptions():
                    assert self.K50 > 0, "K50 has to be positive."
            elif isinstance(self.K50, list):
                with log_exceptions():
                    assert all([True if val > 0 else False for val in self.K50]), \
                        "Subspecies K50 values have to be positive."
            else:  # isinstance(self.K50, tuple)
                with log_exceptions():
                    assert self.K50[0] > 0 and self.K50[1] > 0, \
                        "Mean and std of K50 have to be positive."

            if self.order == 0:
                with log_exceptions():
                    assert not self.is_heterogeneous_K50, \
                        "Heterogeneity in parameter K50 is not allowed for a 0th order process."

            with log_exceptions():
                assert self.nH > 0, "nH has to be positive."

    @classmethod
    def from_string(cls,
                    proc_str: str,
                    /,
                    k: float | int | list[float, ...] | tuple[float, float],
                    *,
                    regulating_species: str | list[str, ...] = None,
                    alpha: float | int | list[float | int, ...] = 1,
                    K50: float | int | list[float | int, ...] | tuple[float | int, float | int] |
                         list[float | int | list[float | int, ...] | tuple[
                             float | int, float | int]] = None,
                    nH: float | int | list[float | int, ...] = None,
                    volume: float | None = None,
                    sep: str = '->') -> Self:
        """ Create a regulated process from a string.

        Parameters
        ----------
        proc_str : str
            A string describing the process in standard chemical notation
            (e.g., 'A + B -> C')
        k : float or int or list of floats or 2-tuple of floats
            The *microscopic* rate constant for the given process. It is the *basal*
            rate constant in the case of activation (or the minimum `k` value) 
            and the maximum rate constant in the case of repression. 
            If `k` is a float or int, then the process is homogeneous. 
            If `k` is a list, then the population of the reactants 
            constsists of distinct subspecies or subinteractions 
            depending on the order. If `k` is a 2-tuple,
            then the constant is normally-distributed with a mean and standard
            deviation specified in the tuple's elements. Note that `k` cannot
            be zero for this form of regulation.
        regulating_species : str or list of str
            Name of the regulating species.
        alpha : float or int or list[float or int]
            Parameter denoting the degree of activation/repression.

                - 0 <= alpha < 1: repression
                - alpha = 1: no regulation
                - alpha > 1: activation
                
            alpha is a multiplier: in the case of activation, the maximum 
            rate constant value will be `alpha * k`. 
            In the case of repression, the minimum 
            rate constant value will be `alpha * k`. 
        K50 : float or int or list of floats or 2-tuple of floats or list of each of the previous types
            *Microscopic* constant that corresponds to the number of
            `regulating_species` agents that would produce 
            half-maximal activation/repression. 
            Heterogeneity in this parameter is determined by the type of `K50`,
            using the same rules as for parameter `k`.
        nH : float or int or list[float or int]
            Hill coefficient for the given process. Indicates the degree of 
            cooperativity in the regulatory interaction.
        volume : float, default : None, optional
            The volume *in liters* of the compartment in which the processes
            are taking place.
        sep : str, default: '->'
            Specifies the characters that distinguish the reactants from the
            products. The default is '->'. The code also treats `-->` as a
            default, if it's present in `proc_str`.

        Notes
        -----
        - Species names should not contain spaces, dashes, and
          should start with a non-numeric character.
        - Zeroth order processes should be specified by an empty space or 'None'.

        Examples
        --------
        >>> RegulatedProcess.from_string("A -> X", k=0.2, regulating_species='X', alpha=2, K50=10, nH=1)
        >>> RegulatedProcess.from_string("A -> X", k=0.3, regulating_species='X', alpha=0.5, K50=[10, 15], nH=2)
        >>> RegulatedProcess.from_string("A + B -> X", k=0.5, regulating_species='B, X', alpha=[2, 0], K50=[(15, 5), [10, 15]], nH=[1, 2])
        """
        sep = '-->' if '-->' in proc_str else sep
        if sep not in proc_str:
            logger.error(f"Cannot find separator '{sep}' in process string: '{proc_str}'")
            raise Exception("Cannot distinguish the reactants from the products.\n"
                            "Please use the *sep* keyword: e.g. sep='->'.")

        lhs, rhs = proc_str.strip().split(sep)  # Left- and Right-hand sides of process
        lhs_terms = lhs.split('+')  # Separate the terms on the left-hand side
        rhs_terms = rhs.split('+')  # Separate the terms on the right-hand side

        return cls(reactants=cls._to_dict(lhs_terms),
                   products=cls._to_dict(rhs_terms),
                   k=k,
                   regulating_species=regulating_species,
                   alpha=alpha,
                   K50=K50,
                   nH=nH,
                   volume=volume)

    def __repr__(self):
        repr_k = macro_to_micro(self.k, self.volume, self.order, inverse=True) if self.volume is not None else self.k
        repr_K50 = macro_to_micro(self.K50, self.volume, inverse=True) if self.volume is not None else self.K50
        return f"RegulatedProcess Object: " \
               f"RegulatedProcess.from_string('{self._str.split(',')[0]}', " \
               f"k={repr_k}, " \
               f"regulating_species='{self.regulating_species}', " \
               f"alpha={self.alpha}, " \
               f"K50={repr_K50}, " \
               f"nH={self.nH}, " \
               f"volume={self.volume})"

    def __str__(self):
        if isinstance(self.regulating_species, list):
            K50_het_str = ""
            for i, sp in enumerate(self.regulating_species):
                if isinstance(self.K50[i], (float, int)):
                    K50_het_str += f"Homogeneous process with respect to species {sp} K50. "
                elif isinstance(self.K50[i], list):
                    K50_het_str += f"Heterogeneous process with respect to species {sp} K50 " \
                                   f"with {len(self.K50[i])} distinct subspecies. "
                else:
                    K50_het_str += f"Heterogeneous process with normally-distributed " \
                                   f"species {sp} K50 with mean {self.K50[i][0]} and " \
                                   f"standard deviation {self.K50[i][1]}. "
        else:
            if isinstance(self.K50, (float, int)):
                K50_het_str = "Homogeneous process with respect to K50."
            elif isinstance(self.K50, list):
                K50_het_str = f"Heterogeneous process with respect to K50 " \
                              f"with {len(self.K50)} distinct subspecies."
            else:
                K50_het_str = f"Heterogeneous process with normally-distributed K50 with " \
                              f"mean {self.K50[0]} and standard deviation {self.K50[1]}."

        return super().__str__() + f" Regulating Species: {self.regulating_species}, " \
                                   f"alpha = {self.alpha}, nH = {self.nH}, " \
                                   f"K50 = {self.K50}, {K50_het_str}"

    def __eq__(self, other):
        if isinstance(other, RegulatedProcess):
            is_equal = (self.k == other.k and
                        self.order == other.order and
                        self.reactants == other.reactants and
                        self.products == other.products and
                        self.regulating_species == other.regulating_species and
                        self.alpha == other.alpha and
                        self.K50 == other.K50 and
                        self.nH == other.nH and
                        self.species == other.species and
                        self.volume == other.volume)
            return is_equal
        elif isinstance(other, str):
            return self._str == other or self._str.replace(' ', '') == other
        else:
            logger.info(f"{type(self)} and {type(other)} are instances of different classes.")
            return False

    def __hash__(self):
        return hash(self._str)


class RegulatedMichaelisMentenProcess(RegulatedProcess):
    """ Define a process that is regulated and obeys Michaelis-Menten kinetics.

    This class allows a Michaelis-Menten Process to be defined
    in terms of how it is regulated.
    If there is only one regulating species, then the parameters have the same
    type as would be expected for a homogeneous/heterogeneous process.
    If there are multiple regulating species, then all parameters are a list
    of their expected type, with the length of the list being equal to the
    number of regulating species.

    The class-specific attributes (except for `k`, which requires some
    additional notes) are listed below.

    Attributes
    ----------
    k : float or int or list of floats or 2-tuple of floats
        The *microscopic* rate constant for the given process. It is the *basal*
        rate constant in the case of activation (or the minimum `k` value)
        and the maximum rate constant in the case of repression.
    regulating_species : str or list of str
        Name of the regulating species. Multiple species can be specified as
        comma-separated in a string or a list of strings with the species names.
    alpha : float or int or list[float or int]
        Parameter denoting the degree of activation/repression.

            - 0 <= alpha < 1: repression
            - alpha = 1: no regulation
            - alpha > 1: activation

        alpha is a multiplier: in the case of activation, the maximum
        rate constant value will be `alpha * k`.
        In the case of repression, the minimum
        rate constant value will be `alpha * k`.
    K50 : float or int or list of floats or 2-tuple of floats or list[float or int or list of floats or 2-tuple of floats]
        *Microscopic* constant that corresponds to the number of
        `regulating_species` agents that would produce
        half-maximal activation/repression.
        Heterogeneity in this parameter is determined by the type of `K50`,
        using the same rules as for parameter `k`.
    nH : float or int or list[float or int]
        Hill coefficient for the given process. Indicates the degree of
        cooperativity in the regulatory interaction.
    is_heterogeneous_K50 : bool or list of bool
        Denotes if the parameter `K50` exhibits heterogeneity
        (distinct subspecies/interactions or normally-distributed).
    regulation_type : str or list of str
        The type of regulation for this process based on the value of alpha:
        'activation' or 'repression' or 'no regulation'.
    catalyst : str
        Name of the species acting as a catalyst for this process.
    Km : float or int or list of floats or 2-tuple of floats
        *Microscopic* Michaelis constant. Corresponds to the number
        of `catalyst` agents that would produce half-maximal activity.
        Heterogeneity in this parameter is determined by the type of `K50`,
        using the same rules as for parameter `k`.
    is_heterogeneous_Km : bool
        Denotes if the parameter `Km` exhibits heterogeneity
        (distinct subspecies/interactions or normally-distributed).

    Notes
    -----
    Currently only implemented for 1st order processes. 0th order processes
    cannot obey Michaelis-Menten kinetics and 2nd order Michaelis-Menten
    processes are not implemented yet.
    """

    def __init__(self,
                 /,
                 reactants: dict[str, int],
                 products: dict[str, int],
                 k: float | int | list[float, ...] | tuple[float, float],
                 *,
                 regulating_species: str | list[str, ...],
                 alpha: float | int | list[float | int, ...],
                 K50: float | int | list[float | int, ...] | tuple[float | int, float | int] |
                      list[float | int | list[float | int, ...] | tuple[float | int, float | int]],
                 nH: float | int | list[float | int, ...],
                 catalyst: str,
                 Km: float | int | list[float | int, ...] | tuple[float | int, float | int],
                 volume: float | None = None):

        self.catalyst = catalyst
        self.Km = Km
        self.is_heterogeneous_Km = False if isinstance(self.Km, (int, float)) else True

        super().__init__(reactants=reactants,
                         products=products,
                         k=k,
                         regulating_species=regulating_species,
                         alpha=alpha,
                         K50=K50,
                         nH=nH,
                         volume=volume)

        super()._validate_reg_params()

        assert self.order != 0, "A 0th order process has no substrate for a catalyst " \
                                "to act on, therefore it cannot follow Michaelis-Menten kinetics."
        if self.order == 2:
            logger.error("2nd order Michaelis-Menten processes are not currently supported.")
            raise NotImplementedError

        if self.volume is not None:  # Convert macroscopic to microscopic Km value
            self.Km = macro_to_micro(Km, self.volume)

        self.species.add(self.catalyst)
        self._str += f", catalyst = {self.catalyst}, Km = {self.Km}"

    @classmethod
    def from_string(cls,
                    proc_str: str,
                    /,
                    k: float | int | list[float, ...] | tuple[float, float],
                    *,
                    regulating_species: str | list[str, ...] = None,
                    alpha: float | int | list[float | int, ...] = 1,
                    K50: float | int | list[float | int, ...] | tuple[float | int, float | int] |
                         list[float | int | list[float | int, ...] | tuple[
                             float | int, float | int]] = None,
                    nH: float | int | list[float | int, ...] = None,
                    catalyst: str = None,
                    Km: float | int | list[float | int, ...] | tuple[
                        float | int, float | int] = None,
                    volume: float | None = None,
                    sep: str = '->') -> Self:
        """ Create a regulated Michaelis-Menten process from a string.

        Parameters
        ----------
        proc_str : str
            A string describing the process in standard chemical notation
            (e.g., 'A + B -> C')
        k : float or int or list of floats or 2-tuple of floats
            The *microscopic* rate constant for the given process. It is the *basal*
            rate constant in the case of activation (or the minimum `k` value)
            and the maximum rate constant in the case of repression.
            If `k` is a float or int, then the process is homogeneous.
            If `k` is a list, then the population of the reactants
            constsists of distinct subspecies or subinteractions
            depending on the order. If `k` is a 2-tuple,
            then the constant is normally-distributed with a mean and standard
            deviation specified in the tuple's elements. Note that `k` cannot
            be zero for this form of regulation.
        regulating_species : str or list of str
            Name of the regulating species.
        alpha : float or int or list[float or int]
            Parameter denoting the degree of activation/repression.

            - 0 <= alpha < 1: repression
            - alpha = 1: no regulation
            - alpha > 1: activation

            alpha is a multiplier: in the case of activation, the maximum
            rate constant value will be `alpha * k`.
            In the case of repression, the minimum
            rate constant value will be `alpha * k`.
        K50 : float or int or list of floats or 2-tuple of floats or list of each of the previous types
            *Microscopic* constant that corresponds to the number of
            `regulating_species` agents that would produce
            half-maximal activation/repression.
            Heterogeneity in this parameter is determined by the type of `K50`,
            using the same rules as for parameter `k`.
        nH : float or int or list[float or int]
            Hill coefficient for the given process. Indicates the degree of
            cooperativity in the regulatory interaction.
        catalyst : str
            Name of species acting as a catalyst.
        Km : float or int or list of floats or 2-tuple of floats
            *Microscopic* Michaelis constant for the process.
            Heterogeneity in this parameter is determined by the type of `Km`,
            using the same rules as for parameter `k`.
        volume : float, default : None, optional
            The volume *in liters* of the compartment in which the processes
            are taking place.
        sep : str, default: '->'
            Specifies the characters that distinguish the reactants from the
            products. The default is '->'. The code also treats `-->` as a
            default, if it's present in `proc_str`.

        Notes
        -----
        - Species names should not contain spaces, dashes, and
          should start with a non-numeric character.
        - Zeroth order processes should be specified by an empty space or 'None'.

        Examples
        --------
        >>> RegulatedMichaelisMentenProcess.from_string("A -> X", k=0.2, regulating_species='X', alpha=2, K50=10, nH=1, catalyst='E', Km=15)
        >>> RegulatedMichaelisMentenProcess.from_string("A -> X", k=0.3, regulating_species='A', alpha=0.5, K50=[10, 15], nH=2, catalyst='C', Km=5)
        """
        sep = '-->' if '-->' in proc_str else sep
        if sep not in proc_str:
            logger.error(f"Cannot find separator '{sep}' in process string: '{proc_str}'")
            raise Exception("Cannot distinguish the reactants from the products.\n"
                            "Please use the *sep* keyword: e.g. sep='->'.")

        lhs, rhs = proc_str.strip().split(sep)  # Left- and Right-hand sides of process
        lhs_terms = lhs.split('+')  # Separate the terms on the left-hand side
        rhs_terms = rhs.split('+')  # Separate the terms on the right-hand side

        return cls(reactants=cls._to_dict(lhs_terms),
                   products=cls._to_dict(rhs_terms),
                   k=k,
                   regulating_species=regulating_species,
                   alpha=alpha,
                   K50=K50,
                   nH=nH,
                   catalyst=catalyst,
                   Km=Km,
                   volume=volume)

    def __repr__(self):
        repr_k = macro_to_micro(self.k, self.volume, self.order, inverse=True) if self.volume is not None else self.k
        repr_K50 = macro_to_micro(self.K50, self.volume, inverse=True) if self.volume is not None else self.K50
        repr_Km = macro_to_micro(self.Km, self.volume, inverse=True) if self.volume is not None else self.Km
        return f"RegulatedMichaelisMentenProcess Object: " \
               f"RegulatedMichaelisMentenProcess.from_string('{self._str.split(',')[0]}', " \
               f"k={repr_k}, " \
               f"regulating_species='{self.regulating_species}', " \
               f"alpha={self.alpha}, " \
               f"K50={repr_K50}, " \
               f"nH={self.nH}, " \
               f"catalyst={self.catalyst}, " \
               f"Km={repr_Km}, " \
               f"volume={self.volume})"

    def __str__(self):
        if isinstance(self.regulating_species, list):
            K50_het_str = ""
            for i, sp in enumerate(self.regulating_species):
                if isinstance(self.K50[i], (float, int)):
                    K50_het_str += f"Homogeneous process with respect to species {sp} K50. "
                elif isinstance(self.K50[i], list):
                    K50_het_str += f"Heterogeneous process with respect to species {sp} K50 " \
                                   f"with {len(self.K50[i])} distinct subspecies. "
                else:
                    K50_het_str += f"Heterogeneous process with normally-distributed " \
                                   f"species {sp} K50 with mean {self.K50[i][0]} and " \
                                   f"standard deviation {self.K50[i][1]}. "
        else:
            if isinstance(self.K50, (float, int)):
                K50_het_str = "Homogeneous process with respect to K50."
            elif isinstance(self.K50, list):
                K50_het_str = f"Heterogeneous process with respect to K50 " \
                              f"with {len(self.K50)} distinct subspecies."
            else:
                K50_het_str = f"Heterogeneous process with normally-distributed K50 with " \
                              f"mean {self.K50[0]} and standard deviation {self.K50[1]}."

        if isinstance(self.Km, (float, int)):
            Km_het_str = "Homogeneous process with respect to Km."
        elif isinstance(self.k, list):
            Km_het_str = f"Heterogeneous process with respect to Km " \
                         f"with {len(self.Km)} distinct subspecies."
        else:
            Km_het_str = f"Heterogeneous process with normally-distributed Km with " \
                         f"mean {self.Km[0]} and standard deviation {self.Km[1]}."

        return super().__str__() + f" Regulating Species: {self.regulating_species}, " \
                                   f"alpha = {self.alpha}, nH = {self.nH}, " \
                                   f"K50 = {self.K50}, {K50_het_str}, " \
                                   f"Catalyst: {self.catalyst}, " \
                                   f"Km = {self.Km}, {Km_het_str}"

    def __eq__(self, other):
        if isinstance(other, RegulatedMichaelisMentenProcess):
            is_equal = (self.k == other.k and
                        self.order == other.order and
                        self.reactants == other.reactants and
                        self.products == other.products and
                        self.regulating_species == other.regulating_species and
                        self.alpha == other.alpha and
                        self.K50 == other.K50 and
                        self.nH == other.nH and
                        self.catalyst == other.catalyst and
                        self.Km == other.Km and
                        self.species == other.species and
                        self.volume == other.volume)
            return is_equal
        elif isinstance(other, str):
            return self._str == other or self._str.replace(' ', '') == other
        else:
            logger.info(f"{type(self)} and {type(other)} are instances of different classes.")
            return False

    def __hash__(self):
        return hash(self._str)


class NullSpeciesNameError(Exception):
    """ Error when the species name is an empty string. """

    def __str__(self):
        return "A species name cannot be an empty string."


def update_all_species(procs: tuple[Process, ...]) -> tuple[set[str], dict[str, list[Process]], dict[str, list[Process]]]:
    """ Categorize all species in a list of processes.

    Extract all species from a list of processes. Then categorize each of them
    as a reactant or product and list the process(es) it takes part in.

    Parameters
    ----------
    procs : tuple
        A tuple of objects of type `Process` or its subclasses.

    Returns
    -------
    tuple
        all_species : set of strings
            A set of all species present in the processes.
        procs_by_reactant : dict
            A dictionary whose keys are the species that are
            reactants in one or more processes. The value for each
            key is a list of processes.
        procs_by_product : dict
            A dictionary whose keys are the species that are
            products in one or more processes. The value for each
            key is a list of processes.
    """
    procs_list = list(procs)
    for proc in procs_list:
        # For a reversible process, replace it with separate instances
        # of Process objects representing the forward and reverse reactions.
        if isinstance(proc, ReversibleProcess):
            forward_proc = Process(proc.reactants, proc.products, proc.k)
            reverse_proc = Process(proc.products, proc.reactants, proc.k_rev)
            procs_list.remove(proc)
            procs_list.extend([forward_proc, reverse_proc])

    with log_exceptions():
        assert len(set(procs_list)) == len(procs_list), \
            "WARNING: Duplicate processes found. Examine the list of processes to resolve this."

    all_species, rspecies, pspecies = set(), set(), set()
    procs_by_reactant, procs_by_product = dict(), dict()

    for proc in procs_list:
        if isinstance(proc, (RegulatedProcess, RegulatedMichaelisMentenProcess)):
            # Add regulating species to the set of all species
            all_species = all_species.union(
                proc.species.union(
                    proc.regulating_species if isinstance(proc.regulating_species, list) else [proc.regulating_species]
                )
            )
        else:
            all_species = all_species.union(proc.species)

        rspecies = rspecies.union(proc.reactants)
        pspecies = pspecies.union(proc.products)

    # Make a list containing the processes each reactant species takes part in.
    # This will be used when solving the system ODEs.
    for rspec in rspecies:
        if rspec != '':  # omit reactant species parsed from 0th order processes
            procs_by_reactant[rspec] = [proc for proc in procs_list if rspec in proc.reactants]
            # deleted 1st clause in above `if`: `rspec != '' and`

    # Make a list containing the processes each product species takes part in.
    # This will be used for solving the system ODEs.
    for pspec in pspecies:
        if pspec != '':  # omitting product species parsed from degradation processes
            procs_by_product[pspec] = [proc for proc in procs_list if pspec in proc.products]

    return all_species, procs_by_reactant, procs_by_product
