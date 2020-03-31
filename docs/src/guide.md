Guide / Tutorial
================================================================================

Installation and tests
--------------------------------------------------------------------------------

BrkgaMpIpr can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia-repl
pkg> add BrkgaMpIpr
```

BrkgaMpIpr also provides a thorough unit testing that aims to harden and make
the code ready for production environments. From Pkg REPL, just run

```julia-repl
pkg> test BrkgaMpIpr
```

!!! note
    The tests take about 10 minutes, mainly because the permutation path relink.

Although BrkgaMpIpr should work fine on Julia >= 1.2, some tests can fail. This
issue occurs because BrkgaMpIpr uses the JLD package to save the population and
results. JLD uses the HDF5 package, which produces slightly different binaries
of different Julia versions. Although the tests may fail in those cases,
BrkgaMpIpr is functional for regular usage. In the table below, you can see
the testing fails due to JDL binary incompatibility.

| Julia version   | Windows   | Linux   | Mac OS X   |
| :-------------: | :-------: | :-----: | :--------: |
| 1.2             | ![](https://img.shields.io/badge/-Fail-red)   | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |
| 1.3             | ![](https://img.shields.io/badge/-Pass-green) | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |
| 1.4             | ![](https://img.shields.io/badge/-Pass-green) | ![](https://img.shields.io/badge/-Pass-green)   | ![](https://img.shields.io/badge/-Pass-green) |

!!! warning
    Some timing tests may fail when carried out on virtual machines and
    containers. The reason is that in such environments, the code runs much
    slower than on bare metal, and some control loops take much time to finish
    before the time stop.  Usually, the difference is a few seconds, but it is
    enough to break some tests.

!!! warning
    It is a hard test to test algorithms that use random signals. In
    BrkgaMpIpr, the tests are carefully designed to ensure repeatability. For
    that, we use the Mersenne Twister
    [[1]](https://en.wikipedia.org/wiki/Mersenne_Twister)
    [[2]](http://dx.doi.org/10.1145/272991.272995) as our standard random
    generator number engine, particularly the [version that comes with
    Julia](https://docs.julialang.org/en/v1/stdlib/Random/index.html#Random.MersenneTwister).
    However, it may happen that such engine has slightly different
    implementations across platforms and, therefore, the tests may fail. The
    current version was tested on 64-bit platforms (Mac OS X, GNU/Linux, and
    Windows 10).

TL;DR
--------------------------------------------------------------------------------

The best way to keep it short is to look in the
[`examples`](https://github.com/ceandrade/brkga_mp_ipr_julia/tree/master/examples/tsp) folder
on [the git repo.](https://github.com/ceandrade/brkga_mp_ipr_julia)
From [`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_minimal.jl),
which solves the
[Travelling Salesman Problem (TSP)](https://en.wikipedia.org/wiki/Travelling_salesman_problem).
This is a trimmed copy:

```julia
using BrkgaMpIpr

include("tsp_instance.jl")
include("tsp_decoder.jl")

if length(ARGS) < 4
    println("Usage: julia main_minimal.jl <seed> <config-file> " *
            "<num-generations> <tsp-instance-file>")
    exit(1)
end

seed = parse(Int64, ARGS[1])
configuration_file = ARGS[2]
num_generations = parse(Int64, ARGS[3])
instance_file = ARGS[4]

instance = TSP_Instance(instance_file)

brkga_data, control_params = build_brkga(
    instance, tsp_decode!, MINIMIZE, seed, instance.num_nodes,
    configuration_file
)

initialize!(brkga_data)

evolve!(brkga_data, num_generations)

best_cost = get_best_fitness(brkga_data)
@show best_cost
```

You can identify the following basic steps:

1. Create a data structure inherited from `AbstractInstance` to hold
   your input data. This object is passed to the decoder function (example
   [`tsp_instance.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/tsp_instance.jl));

2. Implement a decoder function. This function translates a chromosome (array
   of numbers in the interval [0,1]) to a solution for your problem. The decoder
   must return the solution value or cost to be used as fitness by BRKGA
   (example [`tsp_decoder.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/tsp_decoder.jl));

3. Load the instance and other relevant data;

4. Use `build_brkga` to create a `BrkgaData` that represents
   the internal state of the BRKGA-MP-IPR algorithm;

5. Use `initialize!` to init the BRKGA state;

6. Call `evolve!` to optimize;

7. Call `get_best_fitness` and/or `get_best_chromosome` to
   retrieve the best solution.

[`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_minimal.jl)
provides a very minimal example to understand the necessary steps to use the
BRKGA-MP-IPR framework. However,
[`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/master/examples/tsp/main_complete.jl)
provides a full-featured code, handy for scientific use, such as
experimentation and paper writing. This code allows fine-grained control of
the optimization, shows several features of BRKGA-MP-IPR such as the resets,
chromosome injection, and others. It also logs
all optimization steps, _creating outputs easy to be parsed._ **You should use
this code for serious business and experimentation.**

Getting started
--------------------------------------------------------------------------------

BrkgaMpIpr is pretty simple, and you must provide one required data
structure representing the _problem instance_, and one required _decoder_
function to translate chromosomes to solutions.

Before you go further, please take a look at the
[`examples`](https://github.com/ceandrade/brkga_mp_ipr_julia/tree/v1.0/examples)
folder in [the git repo](https://github.com/ceandrade/brkga_mp_ipr_julia).
We will use parts of that code in this guide. There, we solve the classical
[Traveling Salesman
Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem). Given a
set of cities and the distances between them (full weighted undirect graph),
one must find a minimum-cost tour among all cities, such that each city is
visited only once (i.e., find a Hamiltonian cycle of minimum cost). These are
the files:

- [`tsp_instance.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/tsp_instance.jl):
  contains the input data structures and helper functions;

- [`tsp_decoder.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/tsp_decoder.jl):
  contains the decoder function for TSP;

- [`greedy_tour.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/greedy_tour.jl):
  simple heuristic that computes a greedy tour;

- [`config.conf`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/config.conf):
  example of parameter settings;

- [`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_minimal.jl):
  minimal code useful to understand and test the framework.
  **You should start here!** Please take a look on this file before continue
  this tutorial;

- [`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl):
  full-featured code, handy for scientific use, such as
  experimentation and paper writing. This code allows fine-grained control of
  the optimization, shows several features of BRKGA-MP-IPR such as the
  path-relinking calls, resets, chromosome injection, and others. It also logs
  all optimization steps, _creating outputs easy to be parsed._
  **You should use this code for serious business and experimentation;**

- [`instances`](https://github.com/ceandrade/brkga_mp_ipr_julia//tree/v1.0/examples/tsp/instances):
  folder containing some TSP instances for testing.

When you call
[`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_minimal.jl)
or
[`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl)
without arguments,
they show the usage. For example, assuming you are using a terminal:

```bash
$ julia main_minimal.jl
Usage: julia main_minimal.jl <seed> <config-file> <num-generations> <tsp-instance-file>

$ julia main_complete.jl
Usage:
  main_complete.jl -c <config_file> -s <seed> -r <stop_rule> -a <stop_arg> -t <max_time> -i <instance_file> [--no_evolution]
  main_complete.jl (-h | --help)
```

!!! note
    [`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl)
    uses the [DocOpt package](https://github.com/docopt/DocOpt.jl).
    Please, install it before run this script.

So, this is a possible output whe calling
[`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_minimal.jl):

```bash
$ julia main_minimal.jl 27000001 config.conf 100 instances/brazil58.dat
Reading data...
Building BRKGA data and initializing...
Evolving 100 generations...
best_cost = 37552.0
```

For [`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl),
the output is more verbose, since we want to capture
as much information as possible to do some statistical analysis. The output
should be something close to this:

```bash
$ julia main_complete.jl -c config.conf -s 2700001 -r Generations -a 100 -t 60 -i instances/brazil58.dat
------------------------------------------------------
> Experiment started at 2019-02-13T18:40:11.789
> Instance: instances/brazil58.dat
> Configuration: config.conf
> Algorithm Parameters:
>  - population_size 2000
>  - elite_percentage 0.3
>  - mutants_percentage 0.15
>  - num_elite_parents 2
>  - total_parents 3
>  - bias_type LOGINVERSE
>  - num_independent_populations 3
>  - pr_number_pairs 0
>  - pr_minimum_distance 0.15
>  - pr_type PERMUTATION
>  - pr_selection BESTSOLUTION
>  - alpha_block_size 1.0
>  - pr_percentage 1.0
>  - exchange_interval 200
>  - num_exchange_indivuduals 2
>  - reset_interval 600
> Seed: 2700001
> Stop rule: GENERATIONS
> Stop argument: 100
> Maximum time (s): 60.0
> Number of parallel threads for decoding: 1
------------------------------------------------------

[18:40:11.87] Reading TSP data...
Number of nodes: 58

[18:40:11.906] Generating initial tour...
Initial cost: 30774.0

[18:40:11.909] Building BRKGA data...
New population size: 580

[18:40:12.092] Initializing BRKGA data...

[18:40:12.247] Warming up...

[18:40:12.771] Evolving...
* Iteration | Cost | CurrentTime
* 1 | 30774 | 0.03
* 34 | 30751 | 0.83
* 35 | 30507 | 0.85
* 36 | 30088 | 0.87
* 38 | 30023 | 0.93
* 39 | 29882 | 0.95
* 40 | 29665 | 0.97
* 41 | 29131 | 1.00
* 57 | 28221 | 1.38
* 66 | 28211 | 1.59
* 83 | 28200 | 2.01
* 86 | 28129 | 2.08
* 91 | 28118 | 2.19
[18:40:15.171] End of optimization

Total number of iterations: 100
Last update iteration: 91
Total optimization time: 2.40
Last update time: 2.19
Large number of iterations between improvements: 33
Total path relink time: 0.00
Total path relink calls: 0
Number of homogenities: 0
Improvements in the elite set: 0
Best individual improvements: 0

% Best tour cost: 28118
% Best tour: 22 8 1 30 13 40 25 9 32 20 53 50 4 18 44 24 58 5 27 43 12 57 23 54 55 2 41 35 10 52 51 47 49 3 48 39 29 36 17 26 19 6 28 14 37 34 56 46 15 45 33 21 11 16 42 38 31 7

Instance,Seed,NumNodes,TotalIterations,TotalTime,TotalPRTime,PRCalls,NumHomogenities,NumPRImprovElite,NumPrImprovBest,LargeOffset,LastUpdateIteration,LastUpdateTime,Cost
brazil58.dat,2700001,58,100,2.40,0.00,0,0,0,0,33,91,2.19,28118
```

I hope by now you got your system set up and running. Let's see the essential
details on how to use the BrkgaMpIpr.

First things first: basic data structures and decoder function
--------------------------------------------------------------------------------

All problem information must be encapsulated in a struct inherited from
[`AbstractInstance`](@ref). [`AbstractInstance`](@ref) is an empty-abstract
struct required in the signature of the _decoder_ function (described
further). For example, assume we want to solve the Traveling Salesman
Problem. One possible instance struct could be:

```julia
struct TSP_Instance <: AbstractInstance
    num_nodes::Int64
    distances::Array{Float64}
end
```

So, note that we have the number of nodes/cities`num_nodes`, and the
distance matrix `distances`. We may need some additional code to load the
instance and to deal with the distance matrix (omitted here, see
[`tsp_instance.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/tsp_instance.jl)).

The second and most important requirement is the **decoder function**.
The mandatory signature of the decoder is

```julia
decode!(chromosome::Array{Float64, 1},
        problem_instance::AbstractInstance,
        rewrite::Bool)::Float64
```

`chromosome` is a vector of numbers in the interval [0, 1] to be decoded.
`problem_instance` is the data structure containing information about the
problem. Such data is used by the decoder to build a solution. `rewrite` is
an argument that indicates if the decoder should rewrite the
chromosome, in case of local search / local improvements be performed during the
decoder process. This flag is critical if you intend to use the Implicit Path
Relink (details on [`path_relink!`](@ref)). The decoder must return a
`Float64` that is used as the **fitness** to rank the chromosomes. In
general, fitness is the cost/value of the solution, but you may want to use
it to penalize solutions that violate the problem constraints, for example.

In our TSP example, we have a very simple decoder that generates a permutation
of nodes, and compute the cost of the cycle from that permutation. Note the
used of function `distance` that returns the distance between two nodes and
it is defined on
[`tsp_instance.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/tsp_instance.jl).
Also note that, although we do not rewrite the chromosome in this example,
we need `rewrite` in the function signature.

```julia
function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,
                     rewrite::Bool)::Float64

    permutation = Array{Tuple{Float64, Int64}}(undef, instance.num_nodes)
    for (index, key) in enumerate(chromosome)
        permutation[index] = (key, index)
    end

    sort!(permutation)

    cost = distance(instance, permutation[1][2], permutation[end][2])
    for i in 1:(instance.num_nodes - 1)
        cost += distance(instance, permutation[i][2], permutation[i + 1][2])
    end

    return cost
end
```

With the instance data and the decoder ready, we can build the BRKGA data
structures and perform the optimization.

Building BRKGA-MP-IPR data structures
--------------------------------------------------------------------------------

BrkgaMpIpr framework revolves over a single data structure called
[`BrkgaData`](@ref) that represents the internal state of the BRKGA-MP-IPR
algorithm. Since this structure has no constructor, you must build it using
one of the [`Building functions`](@ref building_funcs). There are two
[`build_brkga`](@ref) methods:

* load the parameters from a file:

```julia
function build_brkga(
    problem_instance::AbstractInstance,
    decode_function!::Function,
    opt_sense::Sense,
    seed::Int64,
    chromosome_size::Int64,
    configuration_file::String,
    evolutionary_mechanism_on::Bool = true
)::Tuple{BrkgaData, ExternalControlParams}
```

* load the parameters from a hand-made parameter object:

```julia
function build_brkga(
    problem_instance::AbstractInstance,
    decode_function!::Function,
    opt_sense::Sense,
    seed::Int64,
    chromosome_size::Int64,
    brkga_params::BrkgaParams,
    evolutionary_mechanism_on::Bool = true,
)::BrkgaData
```

Both methods require a `problem_instance` to be used in the
`decode_function!` as commented before.

!!! note
    To date, there is not an easy and clean way to force a function type, as we
    can do in C and C++, and then, `decode_function!` is declared as a simple,
    high-level `Function` type. However, `decode_function!` must have the
    expected signature as explained before.

You also must indicate whether you are minimizing or maximizing through
optimization [`Sense`](@ref).

A good seed also must be provided for the (pseudo) random number generator
(according to [this paper](http://doi.acm.org/10.1145/1276927.1276928)).
BrkgaMpIpr uses the Mersenne Twister
[[1]](http://dx.doi.org/10.1145/272991.272995)
[[2]](https://en.wikipedia.org/wiki/Mersenne_Twister).

The `chromosome_size` also must be given. It indicates the length of each
chromosome in the population. In general, this size depends on the instance
and how the decoder works.

Another common argument is `evolutionary_mechanism_on` which is enabled by
default. When disabled, no evolution is performed. The algorithm only decodes
the chromosomes and ranks them. On each generation, all population is replaced
excluding the best chromosome. This flag helps on implementations of simple
multi-start algorithms.

As said before, the difference between the two methods is from where the
algorithm's hyper-parameters come from. In the first version, the algorithm
reads the BRKGA, IPR, and extra parameters from a simple text file that looks
like this (see
[`config.conf`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/config.conf)
for detailed example):

```txt
population_size 2000
elite_percentage 0.30
mutants_percentage 0.15
num_elite_parents 2
total_parents 3
bias_type LOGINVERSE
num_independent_populations 3
pr_number_pairs 0
pr_minimum_distance 0.15
pr_type PERMUTATION
pr_selection BESTSOLUTION
alpha_block_size 1.0
pr_percentage 1.0
exchange_interval 200
num_exchange_indivuduals 2
reset_interval 600
```

When reading such file, the algorithm ignores all blank lines, and lines
starting with **#**. During the building process, the building method creates a
[`BrkgaParams`](@ref) object and a [`ExternalControlParams`](@ref) object.
[`BrkgaParams`](@ref) contains all hyper-parameters regarding BRKGA and IPR
methods and is stored in the brand-new [`BrkgaData`](@ref) returned by the
method. [`ExternalControlParams`](@ref) are parameters that can be used
outside the BRKGA-MP-IPR to control several aspects of the optimization. For
instance, interval to apply path relink, reset the population, perform
population migration, among others. Although their presence is required on
the config file, they are not mandatory to the BRKGA-MP-IPR itself.

In the second method, we assume we already have a [`BrkgaParams`](@ref)
object, and we just pass it directly to the function. Note that such param
data is deep-copied inside [`BrkgaData`](@ref).

Let's take a look in the example from
[`main_minimal.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_minimal.jl):

```julia
seed = parse(Int64, ARGS[1])
configuration_file = ARGS[2]
num_generations = parse(Int64, ARGS[3])
instance_file = ARGS[4]

instance = TSP_Instance(instance_file)

brkga_data, control_params = build_brkga(
    instance, tsp_decode!, MINIMIZE, seed, instance.num_nodes,
    configuration_file
)
```

This code gets some arguments from the command line and loads a TSP instance.
After that, it builds `brkga_data`. Note that we pass the `instance` data and
the `tsp_decode!` function. Since we are looking for cycles of minimum cost,
we ask for the algorithm `MINIMIZE`. The starting seed is also given. Since
`tsp_decode!` considers each chromosome key as a node/city, the length of the
chromosome must be the number of nodes, i.e., `instance.num_nodes`. Finally,
we also pass the configuration file.

Let's take a look in a more elaborated example
([`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl)):

```julia
brkga_params, control_params = load_configuration(configuration_file)
...
brkga_params.population_size = min(brkga_params.population_size,
                                   10 * instance.num_nodes)
...
brkga_data = build_brkga(instance, tsp_decode!, MINIMIZE, seed,
                         instance.num_nodes, brkga_params, perform_evolution)
```

Here, we first load the configuration file using the helper function
[`load_configuration`](@ref). Then, we modify the population size to be the
minimum between the original parameter and 10x the number of nodes (_making
the population size proportional to the instance size is a very common
strategy used in genetic algorithms_). Then, we call [`build_brkga`](@ref)
using the param data instead of the configuration file (there is an
additional parameter `perform_evolution` to tune the evolution either on or
off, not necessary at this moment).

Now, we have a [`BrkgaData`](@ref) which will be used in all other functions
during the optimization. Note that we can build several [`BrkgaData`](@ref)
objects using different parameters, decoders, or instance data. These
structures can be evolved in parallel and mixed-and-matched at your will.
Each one holds a self-contained BRKGA state including populations, fitness
information, and a state of the random number generator.

Initialization and Warm-start solutions
--------------------------------------------------------------------------------

Before starting the optimization, we need to initialize [`BrkgaData`](@ref)
using [`initialize!`](@ref) function. This procedure initializes the
populations and others data structures of the BRKGA. If an initial population
(warm start) is supplied, the initialization method completes the remaining
individuals, if they do not exist. This method also performs the initial
decoding of the chromosomes. Therefore, depending on the decoder
implementation, this can take a while, and you may want to time such
procedure. The syntax is pretty straightforward:

```julia
initialize!(brkga_data)
```

!!! warning
    `initialize!` must be called before any optimization methods.

!!! warning
    BrkgaMpIpr performs the decoding of each chromosome in parallel if
    multi-thread is enabled. Therefore, **we must guarantee that the decoder is
    THREAD-SAFE.** If such property cannot be held, we suggest using a single
    thread by setting the environmental variable `JULIA_NUM_THREADS = 1`
    [(see Julia Parallel Computing)]
    (https://docs.julialang.org/en/v1/manual/parallel-computing).

### Warm-start solutions

One good strategy is to bootstrap the main optimization algorithm with good
solutions from fast heuristics
[[1](http://dx.doi.org/10.1002/net.21685),
[2](http://dx.doi.org/10.1016/j.ejor.2017.10.045),
[3](http://dx.doi.org/10.1016/j.ejor.2017.10.045)]
or even from relaxations of integer linear programming models
[[4]](http://dx.doi.org/10.1162/EVCO_a_00138).

To do it, you must set these initial solutions before call
[`initialize!`](@ref). Since BRKGA-MP-IPR does not know the problem
structure, you must _encode_ the warm-start solution as chromosomes (vectors
in the interval [0, 1]). In other words, you must do the inverse process that
`decode!` does. For instance, this is a piece of code from
[`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl)
showing this process:

```julia
initial_cost, initial_tour = greedy_tour(instance)
...
keys = sort(rand(instance.num_nodes))
initial_chromosome = zeros(instance.num_nodes)
for i in 1:instance.num_nodes
    initial_chromosome[initial_tour[i]] = keys[i]
end
...
set_initial_population!(brkga_data, [initial_chromosome])
initialize!(brkga_data)
```

Here, we create one incumbent solution using the greedy heuristic
`greedy_tour()` found in
[`greedy_tour.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/greedy_tour.jl).
It gives us `initial_tour` which is a sequence of nodes to be visited. In the
next four lines, we encode `initial_tour`. First, we create a vector of
sorted random `keys`. Note that this is the same order that `tsp_decode!`
uses. We then create the `initial_chromosome`, and fill it up with `keys`
according to the nodes' order in `initial_tour`. Finally, we use
[`set_initial_population!`](@ref) to assign the incumbent to the initial
population. Note that `initial_chromosome` in between braces because
[`set_initial_population!`](@ref) takes a vector of chromosomes. See its
signature:

```julia
set_initial_population!(brkga_data::BrkgaData,
                        chromosomes::Array{Array{Float64, 1}, 1})
```

Indeed, you can have as much warm-start solutions as you like, limited to the
size of the population. Just remember:

!!! warning
    `set_initial_population!` must be called **BEFORE** `initialize!`.

Optimization time: evolving the population
--------------------------------------------------------------------------------

Once all data is set up, it is time to evolve the population and perform
other operations like path-relinking, shaking, migration, and others. The
call is pretty simple:

```julia
evolve!(brkga_data::BrkgaData, num_generations::Int64 = 1)
```

[`evolve!`](@ref) evolves all populations for `num_generations`.

For example, in `main_minimal.jl`, we just evolve the population for a given
number of generations directly and then extract the best solution cost.

```julia
evolve!(brkga_data, num_generations)
best_cost = get_best_fitness(brkga_data)
```

On [`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl),
we have fine-grained control on the optimization. There, we have a main loop
that evolves the population one generation at a time and performs several
operations as to hold the best solution, to check whether it is time for path
relink, population reset, among others. The advantage of that code is that we
can track all optimization details.

!!! warning
    Again, the decoding of each chromosome is done in parallel if
    multi-thread is enabled. Therefore, **we must guarantee that the decoder is
    THREAD-SAFE.** If such property cannot be held, we suggest using a single
    thread by setting the environmental variable `JULIA_NUM_THREADS = 1`
    [(see Julia Parallel Computing)]
    (https://docs.julialang.org/en/v1/manual/parallel-computing).

Accessing solutions/chromosomes
--------------------------------------------------------------------------------

Since Julia does not offer encapsulation mechanisms to keep data private
within data structures, you can access all chromosomes, fitness, and other
data members directly from [`BrkgaData`](@ref). **However, we do not recommend
that, unless you are sure what you are doing.** So, BrkgaMpIpr offers some
helper functions.

Usually, we want to access the best chromosome after some iterations. You can
use the companion functions:

```julia
get_best_fitness(brkga_data::BrkgaData)::Float64
```

```julia
get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}
```

[`get_best_fitness`](@ref) returns the value/fitness of the best chromosome
across all populations.

[`get_best_chromosome`](@ref) returns a _copy_ of the best chromosome across
all populations. You may want to extract an actual solution from such
chromosome, i.e., to apply a decoding function that returns the actual
solution instead only its value.

You may also want to get a copy of specific chromosome for a given population
using [`get_chromosome`](@ref).

```julia
get_chromosome(brkga_data::BrkgaData,
               population_index::Int64,
               position::Int64)::Array{Float64, 1}
```

For example, you can get the 3rd best chromosome from the 2nd population using

```julia
third_best = get_chromosome(brkga_data, 2, 3)
```

Now, suppose you get such chromosome or chromosomes and apply a quick local
search procedure on them. It may be useful to reinsert such new solutions in
the BRKGA population for the next evolutionary cycles. You can do that using
[`inject_chromosome!`](@ref).

```julia
inject_chromosome!(brkga_data::BrkgaData,
                   chromosome::Array{Float64, 1},
                   population_index::Int64,
                   position::Int64,
                   fitness::Float64 = Inf)
```

Note that the chromosome is put in a specific position of a given population.
If you do not provide the fitness, [`inject_chromosome!`](@ref) will decode the
injected chromosome. For example, the following code injects a random chromosome
`keys` into the population #1 in the last position (`population_size`),
i.e., it will replace the worst solution:

```julia
keys = sort(rand(instance.num_nodes))
inject_chromosome!(brkga_data, keys, 1, brkga_data.params.population_size)
```

Implicit Path Relink
--------------------------------------------------------------------------------

The Implicit Path Relinking (IPR) is a nice addition to the standard BRKGA
framework, and it provides an excellent way to create hybrid heuristics and
push the optimization further. The good thing about IPR is that you do not
need to worry about the path relink implementation, which can be long and
tedious if done by hand or customized per problem.

BrkgaMpIpr provides a friendly interface to use IPR directly from the
BRKGA population, and you only must provide a few functions and arguments to
have a Path Relink algorithm ready to go. This is the main signature of
[`path_relink!`](@ref)

```julia
path_relink!(brkga_data::BrkgaData,
             pr_type::PathRelinkingType,
             pr_selection::PathRelinkingSelection,
             compute_distance::Function,
             affect_solution::Function,
             number_pairs::Int64,
             minimum_distance::Float64,
             block_size::Int64 = 1,
             max_time::Float64 = 0.0,
             percentage::Float64 = 1.0
)::PathRelinkingResult
```

The first argument is the [`BrkgaData`](@ref) as usual. The 2nd argument
defines the type of implicit path relink to be performed
([`PathRelinkingType`](@ref)). The `DIRECT` path relink exchanges the keys of
two chromosomes directly, and it is usually more suitable to or threshold
representations, i.e., where the key values are used to some kind of
discretization, such as " if x < 0.5, then 0, otherwise 1." The `PERMUTATION`
path relink switches the order of a key according to its position in the
other chromosome. Usually, this kind of path relink is more suitable to
permutation representations, where the chromosome induces an order or
permutation. For example, chromosome `[0.4, 0.7, 0.1]` may induce the
increasing order `(3, 1, 2)`. More details about threshold and permutation
representations in [this paper](http://dx.doi.org/xxx).

[`PathRelinkingSelection`](@ref) defines how the algorithm picks the
chromosomes for relinking. `BESTSOLUTION` selects, in the order, the best
solution of each population. `RANDOMELITE` chooses uniformly random solutions
from the elite sets.

The next argument is a function to compute the distance between two
chromosomes such signature must be

```julia
compute_distance(vector1::Array{Float64, 1},
                 vector2::Array{Float64, 1})::Float64
```

If the value returned by `compute_distance()` is greater than or equal to
`minimum_distance`, the algorithm will perform the path relink between the
two chromosomes. Otherwise, it will look for another pair of chromosomes.
The algorithm will try `number_pairs` chromosomes before gives up.
In the presence of multiple populations, the path relinking is performed
between elite chromosomes from different populations, in a circular fashion.
For example, suppose we have 3 populations. The framework performs 3 path
relinkings: the first between individuals from populations 1 and 2, the
second between populations 2 and 3, and the third between populations 3 and 1.
In the case of just one population, both base and guiding individuals are
sampled from the elite set of that population.

Note that in traditional path relink algorithms, `compute_distance()` depends
on the problem structure. On IPR, you can use a generic distance function, or
provide one that incorporates more knowledge about the problem. BrkgaMpIpr
provides a function to compute the (modified) [Hamming
distance](https://en.wikipedia.org/wiki/Hamming_distance) for threshold
representations ([`hamming_distance`](@ref)), and a function that computes
the [Kendall Tau
distance](https://en.wikipedia.org/wiki/Kendall_tau_distance) distance for
permutation representations ([`kendall_tau_distance`](@ref)). Again,
details about threshold and permutation representations in
[this paper](http://dx.doi.org/xxx).

As a simple example, suppose you are using a threshold representation where
each chromosome key can represent one of 3 different values (a ternary
threshold representation). So, one possible way to compute the distance
between two chromosomes can be:

```julia
function value(key::Float64)::Float64
    return key < 0.33 ? 0.0 : (key < 0.66 ? 1.0 : 2.0)
end

function compute_distance(vector1::Array{Float64, 1},
                          vector2::Array{Float64, 1})::Float64
    total = 0.0
    for i in 1:length(vector1)
        total += abs(value(vector1[i]) - value(vector2[i]))
    end
    return total
end
```

To avoid changes that do not lead to new solutions, we must verify if such
key exchanges affect the solution. For that, we must pass a function with the
signature:

```julia
affect_solution(block1::SubArray{Float64, 1},
                block2::SubArray{Float64, 1})::Bool
```

`affect_solution` two gets partial chromosomes/block of genes `block1` and
`block2` and checks whether changing the keys from `block1` to `block2`
affects the solution. For instance, suppose that the alleles/keys are used as
threshold such that values > 0.5 activate a feature. Suppose we have
`block1 = [0.3, 0.4, 0.1]` and `block2 = [0.4, 0.1, 0.2]`. Since all values are
below 0.5, changing the keys from `block1` to `block2` do not change the
solution, and therefore, we can drop such change (and subsequently decoding).
The blocks can hold only one key/allele, sequential key blocks, or even the
whole chromosome. Note that `affect_solution` is crucial to the IPR performance
since this function helps to avoid exploring regions already surveyed. Also,
note that `affect_solution` can incorporate some problem knowledge.

!!! warning
    The current implementation of permutation path relink does not make use
    of `affect_solution`. However, [`path_relink!`](@ref) requires the
    function. You can use the simple lambda function for this one:

    ```julia
    (x, y) -> true
    ```

`block_size` defines the number of keys / size of the chromosome block to be
exchanged during the direct path relink. This parameter is also critical for
IPR performance since it avoids too many exchanges during the path building.
Usually, we can compute this number based on the size of the chromosome by
some factor (`alpha_block_size` in the configuration file), chosen by you.
Again, details here.

!!! note
    Experiments have shown that a good choice is

    ``
    block\_size = alpha\_block\_size \times \sqrt{size~of~chromosome}
    ``

The last two parameters are stopping criteria. The algorithm stops either
when `max_time` seconds is reached or `percentage`% of the path is built.

!!! warning
    IPR is a very time-intensive process. You must set the stopping criteria
    accordingly.

Let's see the example on
[`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl).
Remember, since we are solving the TSP, we want to use the permutation-based
IPR, and the Kendall Tau distance functions.

```julia
result = path_relink!(
    brkga_data,
    brkga_params.pr_type,
    brkga_params.pr_selection,
    kendall_tau_distance,
    affect_solution_kendall_tau,
    brkga_params.pr_number_pairs,
    brkga_params.pr_minimum_distance,
    1, #block_size doesn't not matter for permutation.
    maximum_time - (time() - start_time),
    brkga_params.pr_percentage
)
```

Note that most parameters come from [`BrkgaParams`](@ref). The maximum IPR
time is set to the remaining time for optimization (global `maximum_time`
minus the elapsed time `time() - start_time`.

[`path_relink!`](@ref) returns a [`PathRelinkingResult`](@ref) object which
defines the status of the IPR optimization. These status are described on
[`PathRelinkingResult`](@ref).

!!! note
    The `TOO_HOMOGENEOUS` status is directly linked to the chosen distance
    function and minimum distance. If the minimum distance is too large, IPR
    may not be able to find a pair of chromosomes far enough for path relink.

If the found solution is the best solution found so far, IPR replaces the
worst solution by it. Otherwise, IPR computes the distance between the found
solution and all other solutions in the elite set, and replaces the worst
solution by it if and only if the found solution is, at least,
`minimum_distance` from all them.

### Important notes about IPR

The API will call `decode!()` function always with `writeback = false`. The
reason is that if the decoder rewrites the chromosome, the path between
solutions is lost and inadvertent results may come up. Note that at the end
of the path relinking, the method calls the decoder with `writeback = true`
in the best chromosome found to guarantee that this chromosome is re-written
to reflect the best solution found.

!!! warning
    Make sure your decoder does not rewrite the chromosome when called with
    the argument `writeback = false`.

BrkgaMpIpr [`path_relink!`](@ref) implementation is multi-threaded.
Instead of to build and decode each chromosome one at a time, the method
builds a list of candidates, altering the alleles/keys according to the guide
solution, and then decode all candidates in parallel. Note that
``O(chromosome\_size^2 ~/~ block\_size)`` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

!!! warning
    As it is in [`evolve!()`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1/manual/parallel-computing).

Shaking and Resetting
--------------------------------------------------------------------------------

Sometimes, BRKGA gets stuck, converging to local maxima/minima, for several
iterations. When such a situation happens, it is a good idea to perturb the
population, or even restart from a new one completely new. BrkgaMpIpr
offers [`shake!`](@ref) function, an improved variation of the original
version proposed in [this paper](http://dx.doi.org/10.1016/j.eswa.2019.03.007).

```julia
shake!(brkga_data::BrkgaData,
       intensity::Int64,
       shaking_type::ShakingType,
       population_index::Int64 = Inf64)
```

[`shake!`](@ref) function gets an `intensity` parameter that measures how
many times the perturbation is applied on the elite set for a given
`population_index` (if not given, all populations are shaken). This method
offers two generic/implicit [`ShakingType`](@ref)s.
With [`CHANGE`](@ref ShakingType), direct modifications are done in the
keys/alleles. This kind of shaking is recommended when the chromosome uses
direct or threshold representations. [`SWAP`](@ref ShakingType) exchanges
keys/alleles inducing new permutations. For representational definitions,
please read [this paper](http://dx.doi.org/xxx). For instance, the following
code shakes all populations using 10 swap moves.

```julia
shake!(brkga_data, 10, SWAP)
```

Sometimes, even shaking the populations does not help to escape from local
maxima/minima. So, we need a drastic measure, restarting from scratch the
role population. This can be easily accomplished with [`reset!`](@ref).

```julia
reset!(brkga_data)
```

!!! note
    When using [`reset!`](@ref), all warm-start solutions provided by
    [`set_initial_population!`](@ref) are discarded. You may use
    [`inject_chromosome!`](@ref) to insert those solutions again.

!!! warning
    Again, the decoding of each chromosome is done in parallel if
    multi-thread is enabled. Therefore, **we must guarantee that the decoder is
    THREAD-SAFE.** If such property cannot be held, we suggest using a single
    thread by setting the environmental variable `JULIA_NUM_THREADS = 1`
    [(see Julia Parallel Computing)]
    (https://docs.julialang.org/en/v1/manual/parallel-computing).

Multi-population and migration
--------------------------------------------------------------------------------

Multi-population or _island model_ was introduced in genetic algorithms in
[this paper](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.36.7225).
The idea is to evolve parallel and independent populations and, once a
while, exchange individuals among these populations. In several scenarios,
this approach is very beneficial for optimization.

BrkgaMpIpr is implemented using such island idea from the core. If you
read the guide until here, you may notice that several methods take into
account multiple populations. To use multiple populations, you must set
[`BrkgaParams`](@ref)`.num_independent_populations` with 2 ou more populations,
and build [`BrkgaData`](@ref) from such parameters.

The immigration process is implemented by

```julia
exchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)
```

[`exchange_elite!`](@ref) copies `num_immigrants` from one population to
another, replacing the worst `num_immigrants` individuals from the recipient
population. Note that the migration is done for all pairs of populations.
For instance, the following code exchanges 3 best individuals from
each population:

```julia
exchange_elite!(brkga_data, 3)
```

Simulating the standard BRKGA
--------------------------------------------------------------------------------

Sometimes, it is a good idea to test how the standard BRKGA algorithm
performs for a problem. You can use BrkgaMpIpr framework to quickly
implement and test a standard BRKGA.

First, you must guarantee that, during the crossover, the algorithm chooses
only one elite individual and only one non-elite individual. This is easily
accomplished setting `num_elite_parents = 1` and `total_parents = 2`. Then,
you must set up a bias function that ranks the elite and no-elite individual
according to the original BRKGA bias parameter ``\rho`` (rho).

You can use [`set_bias_custom_function!`](@ref) for that task. The given
function receives the index of the chromosome and returns a ranking for it.
Such ranking is used in the roulette method to choose the individual from
which each allele comes to build the new chromosome. Since we have one two
individuals for crossover in the standard BRKGA, the bias function must
return the probability to one or other individual. In the following code, we
do that with a simple `if...else` lambda function.

```julia
# create brkga_params by hand or reading from a file,
# then set the following by hand.
brkga_params.num_elite_parents = 1
brkga_params.total_parents = 2

rho = 0.75
set_bias_custom_function!(brkga_data, x -> x == 1 ? rho : 1.0 - rho)
initialize!(brkga_data)
```

Here, we first set the `num_elite_parents = 1` and `total_parents = 2` as
explained before. Following, we set a variable `rho = 0.75`. This is the
``\rho`` from standard BRKGA, and you may set it as you wish. Then, we set the
bias function as a very simple lambda function:

```julia
x -> x == 1 ? rho : 1.0 - rho
```

So, if the index of the chromosome is 1 (elite individual), it gets a 0.75
rank/probability. If the index is 2 (non-elite individual), the chromosome
gets 0.25 rank/probability.

!!! note
    All these operations must be done before calling [`initialize!`](@ref).


Reading and writing parameters
--------------------------------------------------------------------------------

Although we can build a [`BrkgaData`](@ref) by set up a BrkgaParams object
manually, the easiest way to do so is to read such parameters from a
configuration file. For this, we can use [`load_configuration`](@ref) that
reads a simple plain text file and returns a tuple of [`BrkgaParams`](@ref)
and [`ExternalControlParams`](@ref) objects. For instance,

```julia
brkga_params, control_params = load_configuration("tuned_conf.txt")
```

The configuration file must be plain text such that contains pairs of
parameter name and value. This file must list all fields from
[`BrkgaParams`](@ref) and [`ExternalControlParams`](@ref), even though you do
not use each one. In
[`examples folder`](https://github.com/ceandrade/brkga_mp_ipr_julia/tree/v1.0/examples)
we have
[`config.conf`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/config.conf)
that looks like this:

```txt
population_size 2000
elite_percentage 0.30
mutants_percentage 0.15
num_elite_parents 2
total_parents 3
bias_type LOGINVERSE
num_independent_populations 3
pr_number_pairs 0
pr_minimum_distance 0.15
pr_type PERMUTATION
pr_selection BESTSOLUTION
alpha_block_size 1.0
pr_percentage 1.0
exchange_interval 200
num_exchange_indivuduals 2
reset_interval 600
```

It does not matter whether we use lower or upper cases. Blank lines and lines
starting with **#** are ignored. The order of the parameters should not
matter either. And, finally, this file should be readble for both C++, Julia,
and Python framework versions.

In some cases, you define some of the parameters at the running time, and you
may want to save them for debug or posterior use. To do so, you can use
[`write_configuration`](@ref), call upon a [`BrkgaParams`](@ref) object or
[`BrkgaData`](@ref) object. For example,

```julia
write_configuration("my_new_parameters.conf", brkga_params, external_params)
# or
write_configuration("crazy_parameters.txt", brkga_data)
```

!!! note
    [`write_configuration`](@ref) rewrites the given file. So, watch out to
    not lose previous configurations.

(Probable Valuable) Tips
--------------------------------------------------------------------------------

### Algorithm warmup

When using Julia code, it is an excellent idea to dry-run all functions you
may use and, mainly, the ones you want to time. The reason is that Julia uses
lazy evaluation when live-compiling the code, i.e., it compiles as it goes.
Another advantage is the memory location effects of our data (principle of
locality), that can be brought closer to the processor (L2/L3 caches) during
the running. Obliviously, this depends on how you implement and use your data
structures.

In [`main_complete.jl`](https://github.com/ceandrade/brkga_mp_ipr_julia/blob/v1.0/examples/tsp/main_complete.jl),
we have the following piece of code to warmup mainly the decoder and other
functions. Note that we just deep-copy `brkga_data`, and then, we may lose
the principle of locality.

```julia
bogus_data = deepcopy(brkga_data)
evolve!(bogus_data, 2)
path_relink!(bogus_data, brkga_params.pr_type, brkga_params.pr_selection,
             (x, y) -> 1.0, (x, y) -> true, 0, 0.5, 1, 10.0, 1.0)
get_best_fitness(brkga_data)
get_best_chromosome(brkga_data)
bogus_data = nothing
```

### Complex decoders and timing

Some problems require complex decoders while for others, the decoder contains
local search procedures, that can be time-consuming. In general, the decoding
is the most time-expensive component of a BRKGA algorithm, and it may skew
some stopping criteria based on running time. Therefore, if your decoder is
time-consuming, it is a good idea to implement a timer or chronometer kind of
thing inside the decoder.

Testing for stopping time uses several CPU cycles, and you need to be careful
when/where to test it, otherwise, you spend all the optimization time doing
system calls to the clock.

IMHO, the most effective way to do it is to test time at the very end of the
decoding. If the current time is larger than the maximum time allowed, simple
return `Inf` or `-Inf` according to your optimization direction. In this way,
we make the solution **invalid** since it violates the maximum time allowed.
The BRKGA framework takes care of the rest.

### Multi-threading

Since [Moore's law](https://en.wikipedia.org/wiki/Moore%27s_law) is not
holding its status anymore, we, simple mortals, must appeal to the wonders of
multi-threading. This paradigm can be tricky to code, and [Amdahl's
law](https://en.wikipedia.org/wiki/Amdahl%27s_law) plays against us.
Several genetic algorithms, and in particular, BRKGA, can use parallel
solution evaluation (or decoding), which makes the use of multi-threading
relatively straightforward. BrkgaMpIpr is not different, and it uses
[Julia multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing)
capabilities to do so.

First, as commented several times in this guide, **the decoder must be
THREAD-SAFE.** So, each thread must have its own read/write data structures
and may share other read-only data. The simplest way to do it is to create
those structures inside the decoder (like most people do). **But be aware**,
this strategy slows down the algorithm significantly depending on the size
and format of the structures, and _I do not recommend it_.

IMHO, the best way to do that is to preallocate the data structure per thread
(using
[`Threads.nthreads()`](https://docs.julialang.org/en/v1/base/multi-threading/#Base.Threads.nthreads)),
and pass them to the decoder through the problem instance. Then, inside the
decoder, you can use
[`Threads.threadid()`](https://docs.julialang.org/en/v1/base/multi-threading/#Base.Threads.threadid)
and recover the memory you want to use.

Let's see a simple example considering the TSP example. `tsp_decode!` uses
a single array to create the permutation of nodes. Let's pre-allocate its
memory per thread. So, in `TSP_Instance`, we pre-allocate copies of such
array, one for each thread:

```julia
using Base.Threads

# Declare as a type to make the code shorter and more readable.
PermutationArray = Array{Tuple{Float64, Int64}, 1}

struct TSP_Instance <: AbstractInstance
    num_nodes::Int64
    distances::Array{Float64}

    # Permutations arrays per thread, to be pre-allocated.
    permutation_per_thread::Array{PermutationArray, 1}

    function TSP_Instance(filename::String)
        #... Code for loading here

        # Allocate the main array to create references for other arrays per thread.
        permutation_per_thread = Array{PermutationArray, 1}(undef, nthreads())

        # Pre-allocate the permutation arrays, one for each thread.
        for i in 1:nthreads()
            permutation_per_thread[i] = PermutationArray(undef, num_nodes)
        end

        new(num_nodes, distances, permutation_per_thread)
    end
end
```

Then, in `tsp_decode!`, we simply refer to the array according to the local
thread ID:

```julia
function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,
                     rewrite::Bool)::Float64
    permutation = instance.permutation_per_thread[threadid()]
    #...
    #...
end
```

Note that to pre-allocate decoding structures inside the object holding the
instance is not the most elegant and decoupled code we can write. However, to
decouple the decoding data from the instance data requires that we pass
another data object to the decoder. To do this explicitly, we may get an
embroidered API. We could do it implicitly, by creating a Singleton object to
hold the decoding data. However, this also reduces (a lot) the clarity and
objectivity of the code. In C++ code, this is much easier accomplished by
creating a Decoder object that can hold data members as much as methods.
Therefore, when creating the Decoder object, we can pre-allocate all data
structures we need.

!!! note
    Pre-allocation and multi-threading only make sense for large data
    structures and time-consuming decoders. Otherwise, the code spends too
    much time on context switching and system calls.

