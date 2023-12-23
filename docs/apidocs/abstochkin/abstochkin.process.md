# {py:mod}`abstochkin.process`

```{py:module} abstochkin.process
```

```{autodoc2-docstring} abstochkin.process
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Process <abstochkin.process.Process>`
  - ```{autodoc2-docstring} abstochkin.process.Process
    :summary:
    ```
* - {py:obj}`ReversibleProcess <abstochkin.process.ReversibleProcess>`
  - ```{autodoc2-docstring} abstochkin.process.ReversibleProcess
    :summary:
    ```
* - {py:obj}`MichaelisMentenProcess <abstochkin.process.MichaelisMentenProcess>`
  - ```{autodoc2-docstring} abstochkin.process.MichaelisMentenProcess
    :summary:
    ```
* - {py:obj}`RegulatedProcess <abstochkin.process.RegulatedProcess>`
  - ```{autodoc2-docstring} abstochkin.process.RegulatedProcess
    :summary:
    ```
* - {py:obj}`RegulatedMichaelisMentenProcess <abstochkin.process.RegulatedMichaelisMentenProcess>`
  - ```{autodoc2-docstring} abstochkin.process.RegulatedMichaelisMentenProcess
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`update_all_species <abstochkin.process.update_all_species>`
  - ```{autodoc2-docstring} abstochkin.process.update_all_species
    :summary:
    ```
````

### API

`````{py:class} Process(reactants: dict, products: dict, k: float | int | list[float, ...] | tuple[float, float], **kwargs)
:canonical: abstochkin.process.Process

```{autodoc2-docstring} abstochkin.process.Process
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.Process.__init__
```

````{py:method} _get_reacts_prods_()
:canonical: abstochkin.process.Process._get_reacts_prods_

```{autodoc2-docstring} abstochkin.process.Process._get_reacts_prods_
```

````

````{py:method} from_string(proc_str: str, /, k: float | int | list[float, ...] | tuple[float, float], *, sep: str = '->', **kwargs) -> typing.Self
:canonical: abstochkin.process.Process.from_string
:classmethod:

```{autodoc2-docstring} abstochkin.process.Process.from_string
```

````

````{py:method} _lsp(kwargs: dict)
:canonical: abstochkin.process.Process._lsp
:staticmethod:

```{autodoc2-docstring} abstochkin.process.Process._lsp
```

````

````{py:method} _to_dict(terms: list) -> dict
:canonical: abstochkin.process.Process._to_dict
:staticmethod:

```{autodoc2-docstring} abstochkin.process.Process._to_dict
```

````

````{py:method} _validate_nums()
:canonical: abstochkin.process.Process._validate_nums

```{autodoc2-docstring} abstochkin.process.Process._validate_nums
```

````

````{py:method} __eq__(other)
:canonical: abstochkin.process.Process.__eq__

````

````{py:method} __hash__()
:canonical: abstochkin.process.Process.__hash__

````

````{py:method} __contains__(item)
:canonical: abstochkin.process.Process.__contains__

```{autodoc2-docstring} abstochkin.process.Process.__contains__
```

````

````{py:method} __repr__()
:canonical: abstochkin.process.Process.__repr__

````

````{py:method} __str__()
:canonical: abstochkin.process.Process.__str__

````

````{py:method} _reconstruct_string()
:canonical: abstochkin.process.Process._reconstruct_string

```{autodoc2-docstring} abstochkin.process.Process._reconstruct_string
```

````

`````

`````{py:class} ReversibleProcess(reactants: dict, products: dict, k: float | int | list[float, ...] | tuple[float, float], k_rev: float | int | list[float, ...] | tuple[float, float])
:canonical: abstochkin.process.ReversibleProcess

Bases: {py:obj}`abstochkin.process.Process`

```{autodoc2-docstring} abstochkin.process.ReversibleProcess
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.ReversibleProcess.__init__
```

````{py:method} from_string(proc_str: str, /, k: float | int | list[float, ...] | tuple[float, float], *, k_rev: float | int | list[float, ...] | tuple[float, float] = 0, sep: str = '<->') -> typing.Self
:canonical: abstochkin.process.ReversibleProcess.from_string
:classmethod:

```{autodoc2-docstring} abstochkin.process.ReversibleProcess.from_string
```

````

````{py:method} __repr__()
:canonical: abstochkin.process.ReversibleProcess.__repr__

````

````{py:method} __str__()
:canonical: abstochkin.process.ReversibleProcess.__str__

````

````{py:method} _reconstruct_string()
:canonical: abstochkin.process.ReversibleProcess._reconstruct_string

```{autodoc2-docstring} abstochkin.process.ReversibleProcess._reconstruct_string
```

````

````{py:method} __eq__(other)
:canonical: abstochkin.process.ReversibleProcess.__eq__

````

````{py:method} __hash__()
:canonical: abstochkin.process.ReversibleProcess.__hash__

````

`````

`````{py:class} MichaelisMentenProcess(reactants: dict, products: dict, k: float | int | list[float | int, ...] | tuple[float | int, float | int], catalyst: str, Km: float | int | list[float | int, ...] | tuple[float | int, float | int])
:canonical: abstochkin.process.MichaelisMentenProcess

Bases: {py:obj}`abstochkin.process.Process`

```{autodoc2-docstring} abstochkin.process.MichaelisMentenProcess
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.MichaelisMentenProcess.__init__
```

````{py:method} from_string(proc_str: str, /, k: float | int | list[float | int, ...] | tuple[float | int, float | int], *, catalyst: str = None, Km: float | int | list[float | int, ...] | tuple[float | int, float | int] = None, sep: str = '->') -> typing.Self
:canonical: abstochkin.process.MichaelisMentenProcess.from_string
:classmethod:

```{autodoc2-docstring} abstochkin.process.MichaelisMentenProcess.from_string
```

````

````{py:method} __str__()
:canonical: abstochkin.process.MichaelisMentenProcess.__str__

````

````{py:method} __repr__()
:canonical: abstochkin.process.MichaelisMentenProcess.__repr__

````

````{py:method} __eq__(other)
:canonical: abstochkin.process.MichaelisMentenProcess.__eq__

````

````{py:method} __hash__()
:canonical: abstochkin.process.MichaelisMentenProcess.__hash__

````

`````

`````{py:class} RegulatedProcess(reactants: dict, products: dict, k: float | int | list[float, ...] | tuple[float, float], regulating_species: str | list[str, ...], alpha: float | int | list[float | int, ...], K50: float | int | list[float | int, ...] | tuple[float | int, float | int] | list[float | int | list[float | int, ...] | tuple[float | int, float | int]], nH: float | int | list[float | int, ...])
:canonical: abstochkin.process.RegulatedProcess

Bases: {py:obj}`abstochkin.process.Process`

```{autodoc2-docstring} abstochkin.process.RegulatedProcess
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.RegulatedProcess.__init__
```

````{py:method} _validate_reg_params()
:canonical: abstochkin.process.RegulatedProcess._validate_reg_params

```{autodoc2-docstring} abstochkin.process.RegulatedProcess._validate_reg_params
```

````

````{py:method} from_string(proc_str: str, /, k: float | int | list[float, ...] | tuple[float, float], *, regulating_species: str | list[str, ...] = None, alpha: float | int | list[float | int, ...] = 1, K50: float | int | list[float | int, ...] | tuple[float | int, float | int] | list[float | int | list[float | int, ...] | tuple[float | int, float | int]] = None, nH: float | int | list[float | int, ...] = None, sep: str = '->') -> typing.Self
:canonical: abstochkin.process.RegulatedProcess.from_string
:classmethod:

```{autodoc2-docstring} abstochkin.process.RegulatedProcess.from_string
```

````

````{py:method} __str__()
:canonical: abstochkin.process.RegulatedProcess.__str__

````

````{py:method} __repr__()
:canonical: abstochkin.process.RegulatedProcess.__repr__

````

````{py:method} __eq__(other)
:canonical: abstochkin.process.RegulatedProcess.__eq__

````

````{py:method} __hash__()
:canonical: abstochkin.process.RegulatedProcess.__hash__

````

`````

`````{py:class} RegulatedMichaelisMentenProcess(reactants: dict, products: dict, k: float | int | list[float, ...] | tuple[float, float], regulating_species: str | list[str, ...], alpha: float | int | list[float | int, ...], K50: float | int | list[float | int, ...] | tuple[float | int, float | int] | list[float | int | list[float | int, ...] | tuple[float | int, float | int]], nH: float | int | list[float | int, ...], catalyst: str, Km: float | int | list[float | int, ...] | tuple[float | int, float | int])
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess

Bases: {py:obj}`abstochkin.process.RegulatedProcess`

```{autodoc2-docstring} abstochkin.process.RegulatedMichaelisMentenProcess
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.RegulatedMichaelisMentenProcess.__init__
```

````{py:method} from_string(proc_str: str, /, k: float | int | list[float, ...] | tuple[float, float], *, regulating_species: str | list[str, ...] = None, alpha: float | int | list[float | int, ...] = 1, K50: float | int | list[float | int, ...] | tuple[float | int, float | int] | list[float | int | list[float | int, ...] | tuple[float | int, float | int]] = None, nH: float | int | list[float | int, ...] = None, catalyst: str = None, Km: float | int | list[float | int, ...] | tuple[float | int, float | int] = None, sep: str = '->') -> typing.Self
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess.from_string
:classmethod:

```{autodoc2-docstring} abstochkin.process.RegulatedMichaelisMentenProcess.from_string
```

````

````{py:method} __str__()
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess.__str__

````

````{py:method} __repr__()
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess.__repr__

````

````{py:method} __eq__(other)
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess.__eq__

````

````{py:method} __hash__()
:canonical: abstochkin.process.RegulatedMichaelisMentenProcess.__hash__

````

`````

`````{py:exception} NullSpeciesNameError()
:canonical: abstochkin.process.NullSpeciesNameError

Bases: {py:obj}`Exception`

```{autodoc2-docstring} abstochkin.process.NullSpeciesNameError
```

```{rubric} Initialization
```

```{autodoc2-docstring} abstochkin.process.NullSpeciesNameError.__init__
```

````{py:method} __str__()
:canonical: abstochkin.process.NullSpeciesNameError.__str__

````

`````

````{py:function} update_all_species(procs: tuple[abstochkin.process.Process, ...]) -> tuple[set, dict, dict]
:canonical: abstochkin.process.update_all_species

```{autodoc2-docstring} abstochkin.process.update_all_species
```
````
