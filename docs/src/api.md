API documentation
================================================================================

Enumerations
--------------------------------------------------------------------------------

```@docs
BiasFunction
PathRelinkingResult
PathRelinkingSelection
PathRelinkingType
Sense
ShakingType
```

Types
--------------------------------------------------------------------------------

```@docs
AbstractInstance
BrkgaData
BrkgaParams
ExternalControlParams
```

I/O functions
--------------------------------------------------------------------------------

```@docs
parse
load_configuration
write_configuration
```

[Building functions](@id building_funcs)
--------------------------------------------------------------------------------

```@docs
build_brkga
initialize!
set_bias_custom_function!
set_initial_population!
```

Population manipulation functions
--------------------------------------------------------------------------------

```@docs
exchange_elite!
inject_chromosome!
reset!
shake!
```

Retrival functions
--------------------------------------------------------------------------------

```@docs
get_best_chromosome
get_best_fitness
get_chromosome
get_current_population
```

Evolution functions
--------------------------------------------------------------------------------

```@docs
evolve!
```

Path relink functions
--------------------------------------------------------------------------------

```@docs
affect_solution_hamming_distance
affect_solution_kendall_tau
hamming_distance
kendall_tau_distance
path_relink!
```

Internals
--------------------------------------------------------------------------------

These types and functions are used internally in the framework. They are not
meant to be used directly.

```@meta
CurrentModule = BrkgaMpIpr
```

### Types

```@docs
DecodeStruct
Population
Triple
```

### Minor helper functions
```@docs
:|
empty_function
find_block_range
swap!
```

### Major helper functions
```@docs
evolve_population!
direct_path_relink!
permutation_based_path_relink!
```

Index
--------------------------------------------------------------------------------

```@index
Pages = ["api.md"]
```