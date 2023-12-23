# {py:mod}`abstochkin.base`

```{py:module} abstochkin.base
```

```{autodoc2-docstring} abstochkin.base
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`AbStochKin <abstochkin.base.AbStochKin>`
  - ```{autodoc2-docstring} abstochkin.base.AbStochKin
    :summary:
    ```
````

### API

`````{py:class} AbStochKin(time_unit='sec')
:canonical: abstochkin.base.AbStochKin

```{autodoc2-docstring} abstochkin.base.AbStochKin
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.base.AbStochKin.__init__
```

````{py:method} add_processes_from_file(filename: str)
:canonical: abstochkin.base.AbStochKin.add_processes_from_file

```{autodoc2-docstring} abstochkin.base.AbStochKin.add_processes_from_file
```

````

````{py:method} extract_process_from_str(process_str)
:canonical: abstochkin.base.AbStochKin.extract_process_from_str

```{autodoc2-docstring} abstochkin.base.AbStochKin.extract_process_from_str
```

````

````{py:method} add_process_from_str(process_str: str, /, k: float | int | list[float | int, ...] | tuple[float | int, float | int], **kwargs)
:canonical: abstochkin.base.AbStochKin.add_process_from_str

```{autodoc2-docstring} abstochkin.base.AbStochKin.add_process_from_str
```

````

````{py:method} add_process(reactants: dict, products: dict, /, k: float | int | list[float | int, ...] | tuple[float | int, float | int], **kwargs)
:canonical: abstochkin.base.AbStochKin.add_process

```{autodoc2-docstring} abstochkin.base.AbStochKin.add_process
```

````

````{py:method} del_process_from_str(process_str: str, /, k: float | int | list[float | int, ...] | tuple[float | int], **kwargs)
:canonical: abstochkin.base.AbStochKin.del_process_from_str

```{autodoc2-docstring} abstochkin.base.AbStochKin.del_process_from_str
```

````

````{py:method} del_process(reactants: dict, products: dict, /, k: float | int | list[float | int, ...] | tuple[float | int, float | int], **kwargs)
:canonical: abstochkin.base.AbStochKin.del_process

```{autodoc2-docstring} abstochkin.base.AbStochKin.del_process
```

````

````{py:method} simulate(/, p0: float | int, t_max: float, dt: int = 0.01, n=100, *, random_seed: int = 19, solve_odes: bool = True, ode_method: str = 'RK45', run: bool = True, show_plots: bool = True, multithreading: bool = True, max_agents_by_species: dict = None, max_agents_multiplier: int = 2, _return_simulation: bool = False)
:canonical: abstochkin.base.AbStochKin.simulate

```{autodoc2-docstring} abstochkin.base.AbStochKin.simulate
```

````

````{py:method} simulate_series_in_parallel(series_kwargs: list[dict, ...], *, max_workers: int = None)
:canonical: abstochkin.base.AbStochKin.simulate_series_in_parallel

```{autodoc2-docstring} abstochkin.base.AbStochKin.simulate_series_in_parallel
```

````

`````
