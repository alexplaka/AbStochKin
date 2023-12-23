# {py:mod}`abstochkin.agentstatedata`

```{py:module} abstochkin.agentstatedata
```

```{autodoc2-docstring} abstochkin.agentstatedata
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`AgentStateData <abstochkin.agentstatedata.AgentStateData>`
  - ```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData
    :summary:
    ```
````

### API

`````{py:class} AgentStateData
:canonical: abstochkin.agentstatedata.AgentStateData

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData
```

````{py:attribute} p_init
:canonical: abstochkin.agentstatedata.AgentStateData.p_init
:type: int
:value: >
   None

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.p_init
```

````

````{py:attribute} max_agents
:canonical: abstochkin.agentstatedata.AgentStateData.max_agents
:type: int
:value: >
   None

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.max_agents
```

````

````{py:attribute} reps
:canonical: abstochkin.agentstatedata.AgentStateData.reps
:type: int
:value: >
   None

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.reps
```

````

````{py:attribute} fill_state
:canonical: abstochkin.agentstatedata.AgentStateData.fill_state
:type: int
:value: >
   None

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.fill_state
```

````

````{py:attribute} asv_ini
:canonical: abstochkin.agentstatedata.AgentStateData.asv_ini
:type: numpy.ndarray
:value: >
   'field(...)'

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.asv_ini
```

````

````{py:attribute} asv
:canonical: abstochkin.agentstatedata.AgentStateData.asv
:type: list[numpy.ndarray]
:value: >
   'field(...)'

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.asv
```

````

````{py:method} __post_init__()
:canonical: abstochkin.agentstatedata.AgentStateData.__post_init__

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.__post_init__
```

````

````{py:method} apply_markov_property(r: int)
:canonical: abstochkin.agentstatedata.AgentStateData.apply_markov_property

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.apply_markov_property
```

````

````{py:method} cleanup_asv()
:canonical: abstochkin.agentstatedata.AgentStateData.cleanup_asv

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.cleanup_asv
```

````

````{py:method} get_vals_o1(r: int, stream: numpy.random.Generator, p_vals: numpy.ndarray, state: int = 1)
:canonical: abstochkin.agentstatedata.AgentStateData.get_vals_o1

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.get_vals_o1
```

````

````{py:method} get_vals_o2(other, r: int, stream: numpy.random.Generator, p_vals: numpy.ndarray, state: int = 1)
:canonical: abstochkin.agentstatedata.AgentStateData.get_vals_o2

```{autodoc2-docstring} abstochkin.agentstatedata.AgentStateData.get_vals_o2
```

````

````{py:method} __str__()
:canonical: abstochkin.agentstatedata.AgentStateData.__str__

````

`````
