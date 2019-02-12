Guide
================================================================================

The BrkgaMpIpr.jl is pretty simple, and you must provide one required data
structure representing the _problem instance_, and one required _decoder_
function to translate chromosomes to solutions.

Before you go further, please take a look at the `examples` folder in
[the git repo](https://github.com/ceandrade/BrkgaMpIpr).
We will use parts of that code in this tutorial. There, we solve
the classical [Traveling Salesman
Problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem). Given a
set of cities and the distances between them (full weighted undirect graph),
one must find a minimum-cost tour among all cities, such that each city is
visited only once (i.e., find a Hamiltonian cycle of minimum cost). These are
the files:

  - `tsp_instance.jl`: contains the input data structures and helper functions;

  - `tsp_decoder.jl`: contains the decoder function for TSP;

  - `greedy_tour.jl`: simple heuristic that computes a greedy tour;

  - `config.conf`: example of parameter settings;

  - `main_minimal.jl`: minimal code useful to understand and test the framework.
    **You should start here!**. Please take a look on this file before continue
    this tutorial;

  - `main_complete.jl`: full-featured code, handy for scientific use, such as
    experimentation and paper writing. This code allows fine-grained control of
    the optimization, show several features of BRKGA-MP-IPR such as the
    path-reliking calls, resets, chromosome injection, and others. It also logs
    all optimization steps, creating outputs easy to be parsed.
    **You should use this code for serious business and experimentation;**

  - `instances`: folder containing some TSP instances for testing.

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
instance and to deal with the distance matrix (omitted here,
see `tsp_instance.jl`).

The second and most important requirement is the **decoder function**.
The mandatory signature of the decoder is

```julia
decode!(chromosome::Array{Float64, 1},
        problem_instance::AbstractInstance,
        rewrite::Bool = true)::Float64
```

`chromosome` is a vector of numbers in the interval [0, 1] to be decoded.
`problem_instance` is the data structure containing information about the
problem. Such data is used by the decoder to build a solution. `rewrite` is
an optional argument that indicates if the decoder should rewrite the
chromosome, in case of local search / local improvements be performed during the
decoder process. This flag is critical if you intend to use the Implicit Path
Relink (details on [`path_relink!`](@ref)). The decoder must return a
`Float64` that is used as the **fitness** to rank the chromosomes. In
general, fitness is the cost/value of the solution, but you may want to use
it to penalize solutions that violate the problem constraints, for example.

In our TSP example, we have a very simple decoder that generates a permutation
of nodes, and compute the cost of the cycle from that permutation (note the
used of function `distance` that returns the distance between two nodes and
it si defined on `tsp_instance.jl`):

```julia
function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,
                     rewrite::Bool = true)::Float64

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

BrkgaMpIpr.jl framework revolves over a single data structure called
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

A good seed also must beprovided for the (pseudo) random number generator
(BrkgaMpIpr.jl uses the Mersenne Twister
[[1]](https://en.wikipedia.org/wiki/Mersenne_Twister)
[[2]](http://dx.doi.org/10.1145/272991.272995)).

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
like this (see `config.conf` for detailed example):

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

Let's take a look in the example from `main_minimal.jl`:

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

Let's take a look in a more elaborated example (`main_complete.jl`):

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
    BrkgaMpIpr.jl performs the decoding of each chromosome in parallel if
    multi-thread is enabled. Therefore, **we must guarantee that the decoder is
    THREAD-SAFE.** If such property cannot be held, we suggest using a single
    thread by setting the environmental variable `JULIA_NUM_THREADS = 1`
    [(see Julia Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

### Warm-start solutions

One good strategy is to bootstrap the main optimization algorithm with good
solutions from fast heuristics
[[1](http://dx.doi.org/10.1002/net.21685),
[2](http://dx.doi.org/10.1016/j.ejor.2017.10.045),
[3](http://dx.doi.org/10.1016/j.ejor.2017.10.045)]
or even from relaxations of integer linear programming models
[[4]](http://dx.doi.org/10.1162/EVCO_a_00138).

To do it, you must set these initial solutions before call [`initialize!`](@ref).
Since BRKGA-MP-IPR does not know the problem structure, you must _encode_ the
warm-start solution as chromosomes (vectors in the interval [0, 1]). In other
words, you must do the inverse process that `decode!` does. For instance,
this is a piece of code from `main_complete.jl` showing this process:

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
`greedy_tour()` (in `greedy_tour.jl`). It gives us `initial_tour` which is a
sequence of nodes to be visited. In the next four lines, we encode
`initial_tour`. First, we create a vector of sorted random `keys`. Note that
this is the same order that `tsp_decode!` uses. We then create the
`initial_chromosome`, and fill it up with `keys` according to the nodes'
order in `initial_tour`. Finally, we use [`set_initial_population!`](@ref) to
assign the incumbent to the initial population. Note that
`initial_chromosome` in between braces because
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

On `main_complete.jl`, we have fine-grained control on the optimization.
There, we have a main loop that evolves the population one generation at a
time and performs several operations as to hold the best solution, to check
whether it is time for path relink, population reset, among others. The
advantage of that code is that we can track all optimization details.

!!! warning
    Again, the decoding of each chromosome is done in parallel if
    multi-thread is enabled. Therefore, **we must guarantee that the decoder is
    THREAD-SAFE.** If such property cannot be held, we suggest using a single
    thread by setting the environmental variable `JULIA_NUM_THREADS = 1`
    [(see Julia Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

Accessing solutions/chromosomes
--------------------------------------------------------------------------------

Since Julia does not offer encapsulation mechanisms to keep data private
within data structures, you can access all chromosomes, fitness, and other
data members directly from [`BrkgaData`](@ref). **However, we do not recommend
that, unless you are sure what you are doing.** So, BrkgaMpIpr.jl offers some
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
`keys` into the population #1 in the last position (`population_size`).

```julia
keys = sort(rand(instance.num_nodes))
inject_chromosome!(brkga_data, keys, 1, brkga_data.params.population_size)
```

Implicit Path Relink
--------------------------------------------------------------------------------


Shaking and Resetting
--------------------------------------------------------------------------------

Sometimes, BRKGA gets stuck, converging to local maxima/minima, for several iterations. When such a situation happens, it is a good idea to perturb the population, or even restart from a new one completely new.

BrkgaMpIpr.jl offers [`shake!`](@ref) function, an improved variation of the original version proposed in
[this paper](http://dx.doi.org/xxx).

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
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

Multi-population and migration
--------------------------------------------------------------------------------

Multi-population or _island model_ was introduced in genetic algorithms in
[this paper](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.36.7225).
The idea is to evolve parallel and independent populations and, once a
while, exchange individuals among these populations. In several scenarios,
this approach is very beneficial for optimization.

BrkgaMpIpr.jl is implemented using such island idea from the core. If you
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
