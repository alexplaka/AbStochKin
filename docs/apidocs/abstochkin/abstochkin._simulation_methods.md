# {py:mod}`abstochkin._simulation_methods`

```{py:module} abstochkin._simulation_methods
```

```{autodoc2-docstring} abstochkin._simulation_methods
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`SimulationMethodsMixin <abstochkin._simulation_methods.SimulationMethodsMixin>`
  - ```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin
    :summary:
    ```
````

### API

`````{py:class} SimulationMethodsMixin
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin
```

````{py:method} _validate_p0()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._validate_p0

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._validate_p0
```

````

````{py:method} _setup_data()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._setup_data

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._setup_data
```

````

````{py:method} _setup_runtime_data()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._setup_runtime_data

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._setup_runtime_data
```

````

````{py:method} _get_fill_state(species_name)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._get_fill_state

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._get_fill_state
```

````

````{py:method} _init_het_metrics()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_het_metrics

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_het_metrics
```

````

````{py:method} _get_o0_trans_p(proc: abstochkin.process.Process | abstochkin.process.RegulatedProcess)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._get_o0_trans_p

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._get_o0_trans_p
```

````

````{py:method} _init_o1_vals(proc: abstochkin.process.Process)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_o1_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_o1_vals
```

````

````{py:method} _init_o1mm_Km_vals(proc: abstochkin.process.MichaelisMentenProcess)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_o1mm_Km_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_o1mm_Km_vals
```

````

````{py:method} _init_o1reg_K50_vals(proc: abstochkin.process.RegulatedProcess)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_o1reg_K50_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_o1reg_K50_vals
```

````

````{py:method} _gen_o1_distinct_subspecies_vals(proc: abstochkin.process.Process | abstochkin.process.MichaelisMentenProcess | abstochkin.process.RegulatedProcess, num_agents: int, *, het_attr: str = 'k', het_attr_idx: int = None)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_o1_distinct_subspecies_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_o1_distinct_subspecies_vals
```

````

````{py:method} _gen_o1_normal_vals(proc: abstochkin.process.Process | abstochkin.process.MichaelisMentenProcess | abstochkin.process.RegulatedProcess, num_agents: int, *, het_attr: str = 'k', het_attr_idx: int = None, output_type: type[float | int] = float, max_tries: int = 1000)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_o1_normal_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_o1_normal_vals
```

````

````{py:method} _init_o2_vals(proc: abstochkin.process.Process)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_o2_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_o2_vals
```

````

````{py:method} _init_o2reg_K50_vals(proc: abstochkin.process.RegulatedProcess)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._init_o2reg_K50_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._init_o2reg_K50_vals
```

````

````{py:method} _gen_o2_distinct_subinteractions_vals(proc: abstochkin.process.Process | abstochkin.process.RegulatedProcess, phm_shape: tuple[int, int], *, het_attr: str = 'k', het_attr_idx: int = None)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_o2_distinct_subinteractions_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_o2_distinct_subinteractions_vals
```

````

````{py:method} _gen_o2_normal_vals(proc: abstochkin.process.Process | abstochkin.process.RegulatedProcess, phm_shape: tuple[int, int], *, het_attr: str = 'k', het_attr_idx: int = None, output_type: type[float | int] = float, max_tries: int = 1000)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_o2_normal_vals

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_o2_normal_vals
```

````

````{py:method} _gen_algo_processes(procs)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_algo_processes

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_algo_processes
```

````

````{py:method} _gen_algo_sequence()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._gen_algo_sequence

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._gen_algo_sequence
```

````

````{py:method} _parallel_run(num_workers=None)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._parallel_run

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._parallel_run
```

````

````{py:method} _sequential_run()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._sequential_run

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._sequential_run
```

````

````{py:method} _repeat_sim(r: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._repeat_sim

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._repeat_sim
```

````

````{py:method} _order_0(proc: abstochkin.process.Process, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_0

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_0
```

````

````{py:method} _order_0_reg(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_0_reg

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_0_reg
```

````

````{py:method} _order_0_reg_gt1(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_0_reg_gt1

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_0_reg_gt1
```

````

````{py:method} _order_1(proc: abstochkin.process.Process, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1
```

````

````{py:method} _order_1_mm(proc: abstochkin.process.MichaelisMentenProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1_mm

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1_mm
```

````

````{py:method} _order_1_reg(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg
```

````

````{py:method} _order_1_reg_gt1(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_gt1

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_gt1
```

````

````{py:method} _order_1_reg_mm(proc: abstochkin.process.RegulatedMichaelisMentenProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_mm

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_mm
```

````

````{py:method} _order_1_reg_mm_gt1(proc: abstochkin.process.RegulatedMichaelisMentenProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_mm_gt1

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_1_reg_mm_gt1
```

````

````{py:method} _order_2(proc: abstochkin.process.Process, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_2

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_2
```

````

````{py:method} _order_2_reg(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_2_reg

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_2_reg
```

````

````{py:method} _order_2_reg_gt1(proc: abstochkin.process.RegulatedProcess, r: int, t: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._order_2_reg_gt1

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._order_2_reg_gt1
```

````

````{py:method} _compute_k_het_metrics(proc: abstochkin.process.Process | abstochkin.process.MichaelisMentenProcess | abstochkin.process.RegulatedProcess | abstochkin.process.RegulatedMichaelisMentenProcess, all_k_vals: numpy.array, t: int, r: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._compute_k_het_metrics

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._compute_k_het_metrics
```

````

````{py:method} _compute_Km_het_metrics(proc: abstochkin.process.MichaelisMentenProcess | abstochkin.process.RegulatedMichaelisMentenProcess, all_Km_vals: numpy.array, t: int, r: int)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._compute_Km_het_metrics

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._compute_Km_het_metrics
```

````

````{py:method} _compute_K50_het_metrics(proc: abstochkin.process.RegulatedProcess | abstochkin.process.RegulatedMichaelisMentenProcess, all_K50_vals: numpy.array, t: int, r: int, *, idx: int = None)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._compute_K50_het_metrics

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._compute_K50_het_metrics
```

````

````{py:method} _o2_get_unique_pairs(i_pairs: numpy.array)
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._o2_get_unique_pairs
:staticmethod:

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._o2_get_unique_pairs
```

````

````{py:method} _compute_trajectory_stats()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._compute_trajectory_stats

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._compute_trajectory_stats
```

````

````{py:method} _compute_het_stats()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._compute_het_stats

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._compute_het_stats
```

````

````{py:method} _post_run_cleanup()
:canonical: abstochkin._simulation_methods.SimulationMethodsMixin._post_run_cleanup

```{autodoc2-docstring} abstochkin._simulation_methods.SimulationMethodsMixin._post_run_cleanup
```

````

`````
