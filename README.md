# AbStochKin: Agent-based Stochastic Kinetics

<p>
  <a href="https://doi.org/10.5281/zenodo.14255157">
    <img src="https://zenodo.org/badge/733779271.svg" alt="DOI">
  </a>
  <a href="https://pypi.org/project/abstochkin/">
    <img src="https://img.shields.io/pypi/v/abstochkin"
         alt="PyPI package version">
  </a>
  <a href="https://www.python.org">
    <img src="https://img.shields.io/pypi/pyversions/abstochkin"
         alt="Python versions">
  </a>
  <a href="https://alexplaka.github.io/AbStochKin">
    <img src="https://img.shields.io/badge/-documentation-blue"
         alt="Documentation">
  </a>
  <a href="https://github.com/alexplaka/abstochkin/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/alexplaka/abstochkin"
         alt="GitHub license">
  </a>
</p>

##### Alternate name: PyStochKin (Particle-based Stochastic Kinetics)

`AbStochKin` is an agent-based (or particle-based) Monte-Carlo simulator of the time
evolution of systems composed of species that participate in coupled processes.
The population of a species is considered as composed of distinct individuals,
termed *agents*, or *particles*.
This allows for the specification of the kinetic parameters describing
the propensity of *each agent* to participate in a given process.

Although the algorithm was originally conceived for simulating biochemical
systems, it is applicable to other disciplines where there is a need to model
how populations change over time and to study the effects of heterogeneity,
or diversity, in the composition of species populations on the dynamics of a system.

## Installation

The `abstochkin` package can be installed via `pip` in an environment with Python 3.12+.

`$ pip install abstochkin`

### Requirements

The package relies only on Python's scientific ecosystem
libraries (`numpy`, `scipy`, `sympy`, `matplotlib`) and
the standard library for implementing the core components of the algorithm.
These requirements can be easily met in any Python (version 3.12+) environment.

## What processes can be modeled?

- Simple processes (0th, 1st, 2nd order).
- Processes obeying Michaelis-Menten kinetics (1st order).
- Processes that are regulated by one or more species through activation or repression (0th, 1st, 2nd order).
- Processes that are regulated *and* obey Michaelis-Menten kinetics (1st order).

## Usage

Here is a simple example of how to run a simulation: consider the process
$A \rightarrow B$, the conversion of agents of species $A$ to agents of species $B$.
Notice that we represent the process in standard chemical notation, therefore
there are 'reactants' and 'products' and each species has a stoichiometric
coefficient associated with it (implied to be $1$ if it is not explicitly written).
The rate constant for this process is specified to be $k=0.2$ and has units of
reciprocal seconds. Here, we assume a homogeneous population; that is, all
agents of species $A$ have the same propensity to 'transition' to species $B$.
Thus, the value $k=0.2$ applies to all $A$ agents when determining the transition
probability within a given time step.

We then run an ensemble of simulations by specifying the initial population
sizes ($A$: $100$ agents, $B$: $0$ agents) and the simulated time of $10$ seconds.
Behind the scenes, default values for unspecified but necessary arguments are used
(specifically, the number of simulations that comprise the ensemble, $n=100$,
and the duration of the fixed time interval for each step in the simulation,
$dt=0.01$ seconds).

```python
from abstochkin import AbStochKin

sim = AbStochKin()
sim.add_process_from_str('A -> B', k=0.2)
sim.simulate(p0={'A': 100, 'B': 0}, t_max=10)
```

When the simulation is completed, the results are presented in graphical form.

## Concurrency

The algorithm performs an ensemble of simulations to obtain the mean time
trajectory of all species and statistical measures of the uncertainty thereof.
To facilitate the rapid execution of the simulation, *multithreading* is enabled
by default. This is done because `numpy`, whose core algorithms can bypass
the Global Interpreter Lock (GIL), is used extensively during the algorithm's
runtime. For instance, the simple usage example presented above uses
multithreading.

When running a series of jobs (each with its own ensemble of simulations)
where a parameter is varied (e.g., a parameter sweep), *process-based
parallellism* can be used. The user does not have to worry about
the details of setting up the code for multiprocessing. Instead, they can simply
call a method of the base class.

```python
from abstochkin import AbStochKin

sim = AbStochKin()
# Define a process that obeys Michaelis-Menten kinetics:
sim.add_process_from_str("A -> B", k=0.3, catalyst='E', Km=10)
# Vary the initial population size of species A:
series_kwargs = [{"p0": {'A': a, 'B': 0, 'E': 10}, "t_max": 10} for a in range(40, 51)]
sim.simulate_series_in_parallel(series_kwargs)
```

## Documentation

See the documentation [here](https://alexplaka.github.io/AbStochKin).

A monograph detailing the theoretical underpinnings of the *Agent-based Kinetics*
algorithm and a multitude of case studies highlighting its use can be found
[here](/docs/Agent-basedKinetics_monograph.pdf).

## Contributing

We welcome any contributions to the project in the form of bug reports,
feature requests, and pull requests. Feel free to contact the core developer
and maintainer at alex dot plaka at alumni dot princeton.edu to introduce
yourself and discuss possible ways to contribute.
