# {py:mod}`abstochkin.simulation`

```{py:module} abstochkin.simulation
```

```{autodoc2-docstring} abstochkin.simulation
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Simulation <abstochkin.simulation.Simulation>`
  - ```{autodoc2-docstring} abstochkin.simulation.Simulation
    :summary:
    ```
````

### API

`````{py:class} Simulation(/, p0: float, t_max: float, dt: int, n: list, processes, *, random_state: int, do_solve_ODEs: bool, ODE_method: str, do_run: bool, show_graphs: bool, use_multithreading: bool, max_agents: dict, max_agents_multiplier: float | int, time_unit: str)
:canonical: abstochkin.simulation.Simulation

Bases: {py:obj}`abstochkin._simulation_methods.SimulationMethodsMixin`

```{autodoc2-docstring} abstochkin.simulation.Simulation
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.simulation.Simulation.__init__
```

````{py:method} run_simulation()
:canonical: abstochkin.simulation.Simulation.run_simulation

```{autodoc2-docstring} abstochkin.simulation.Simulation.run_simulation
```

````

````{py:method} graph_results(/, graphs_to_show=None, species_to_show=None)
:canonical: abstochkin.simulation.Simulation.graph_results

```{autodoc2-docstring} abstochkin.simulation.Simulation.graph_results
```

````

````{py:method} __repr__()
:canonical: abstochkin.simulation.Simulation.__repr__

````

````{py:method} __str__()
:canonical: abstochkin.simulation.Simulation.__str__

````

`````
