# {py:mod}`abstochkin.de_calcs`

```{py:module} abstochkin.de_calcs
```

```{autodoc2-docstring} abstochkin.de_calcs
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`DEcalcs <abstochkin.de_calcs.DEcalcs>`
  - ```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs
    :summary:
    ```
````

### API

`````{py:class} DEcalcs(p0: dict, t_min: int | float, t_max: int | float, processes: list, ode_method: str, time_unit: str)
:canonical: abstochkin.de_calcs.DEcalcs

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs.__init__
```

````{py:method} setup_ODEs(agent_based=True)
:canonical: abstochkin.de_calcs.DEcalcs.setup_ODEs

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs.setup_ODEs
```

````

````{py:method} get_term_multiplier(proc)
:canonical: abstochkin.de_calcs.DEcalcs.get_term_multiplier
:staticmethod:

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs.get_term_multiplier
```

````

````{py:method} solve_ODEs()
:canonical: abstochkin.de_calcs.DEcalcs.solve_ODEs

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs.solve_ODEs
```

````

````{py:method} get_fixed_pts()
:canonical: abstochkin.de_calcs.DEcalcs.get_fixed_pts

```{autodoc2-docstring} abstochkin.de_calcs.DEcalcs.get_fixed_pts
```

````

`````
