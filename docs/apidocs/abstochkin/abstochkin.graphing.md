# {py:mod}`abstochkin.graphing`

```{py:module} abstochkin.graphing
```

```{autodoc2-docstring} abstochkin.graphing
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Graph <abstochkin.graphing.Graph>`
  - ```{autodoc2-docstring} abstochkin.graphing.Graph
    :summary:
    ```
````

### API

`````{py:class} Graph(/, nrows=1, ncols=1, figsize=(5, 5), dpi=300, **kwargs)
:canonical: abstochkin.graphing.Graph

```{autodoc2-docstring} abstochkin.graphing.Graph
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.graphing.Graph.__init__
```

````{py:method} setup_spines_ticks(ax_loc)
:canonical: abstochkin.graphing.Graph.setup_spines_ticks

```{autodoc2-docstring} abstochkin.graphing.Graph.setup_spines_ticks
```

````

````{py:method} plot_ODEs(de_data, *, num_pts: int = 1000, species: list[str] | tuple[str] = (), ax_loc: tuple = ())
:canonical: abstochkin.graphing.Graph.plot_ODEs

```{autodoc2-docstring} abstochkin.graphing.Graph.plot_ODEs
```

````

````{py:method} plot_trajectories(time, data, *, species: list[str] | tuple[str] = (), ax_loc: tuple = ())
:canonical: abstochkin.graphing.Graph.plot_trajectories

```{autodoc2-docstring} abstochkin.graphing.Graph.plot_trajectories
```

````

````{py:method} plot_avg_std(time, data, *, species: list[str] | tuple[str] = (), ax_loc: tuple = ())
:canonical: abstochkin.graphing.Graph.plot_avg_std

```{autodoc2-docstring} abstochkin.graphing.Graph.plot_avg_std
```

````

````{py:method} plot_eta(time, data, *, species: list[str] | tuple[str] = (), ax_loc: tuple = ())
:canonical: abstochkin.graphing.Graph.plot_eta

```{autodoc2-docstring} abstochkin.graphing.Graph.plot_eta
```

````

````{py:method} plot_het_metrics(time, proc_str: tuple[str, str], proc_data: dict, *, het_attr='k', ax_loc: tuple = ())
:canonical: abstochkin.graphing.Graph.plot_het_metrics

```{autodoc2-docstring} abstochkin.graphing.Graph.plot_het_metrics
```

````

````{py:method} savefig(filename, **kwargs)
:canonical: abstochkin.graphing.Graph.savefig

```{autodoc2-docstring} abstochkin.graphing.Graph.savefig
```

````

`````
