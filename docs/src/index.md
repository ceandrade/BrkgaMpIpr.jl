# BrkgaMpIpr.jl Documentation

## Enumerations
```@docs
Sense
BiasFunction
PathRelinkingType
PathRelinkingSelection
PathRelinkingResult
ShakingType
``` 

## Types
```@docs
BrkgaData
AbstractInstance
BrkgaParams
ExternalControlParams
```


## I/O functions
```@doc
parse
load_configuration
write_configuration
```


## Building functions
```@docs
build_brkga
set_bias_custom_function!
set_initial_population!
initialize!
```

## Support functions
```@docs
reset!
exchange_elite!
inject_chromosome!
shake!
```


## Retrival functions
```@docs
get_best_fitness
get_best_chromosome
get_chromosome
get_current_population
```

## Evolution functions
```@docs
evolve!
```

## Path relink functions
```@docs
hamming_distance
affect_solution_hamming_distance
kendall_tau_distance
affect_solution_kendall_tau
path_relink!
```

