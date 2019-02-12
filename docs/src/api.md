API documentation
================================================================================

Enumerations
--------------------------------------------------------------------------------

```@docs
Sense
BiasFunction
PathRelinkingType
PathRelinkingSelection
PathRelinkingResult
ShakingType
```

Types
--------------------------------------------------------------------------------

```@docs
BrkgaData
AbstractInstance
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
set_bias_custom_function!
set_initial_population!
initialize!
```

Support functions
--------------------------------------------------------------------------------

```@docs
reset!
exchange_elite!
inject_chromosome!
shake!
```

Retrival functions
--------------------------------------------------------------------------------

```@docs
get_best_fitness
get_best_chromosome
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
hamming_distance
affect_solution_hamming_distance
kendall_tau_distance
affect_solution_kendall_tau
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
Population
Triple
DecodeStruct
empty_function
```

### Helper functions

```@docs
evolve_population!
direct_path_relink!
permutation_based_path_relink!
swap!
:|
find_block_range
```
