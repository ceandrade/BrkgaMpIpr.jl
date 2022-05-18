################################################################################
# building.jl: Building and parameter setup routines.
#
# (c) Copyright 2022, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: May 18, 2022 by ceandrade
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
################################################################################

"""
    const empty_function() = nothing

Represent an empty function to be used as flag during data and bias function
setups.
"""
const empty_function() = nothing

################################################################################

"""
    build_brkga(problem_instance::AbstractInstance, decode_function!::Function,
                opt_sense::Sense, seed::Int64, chromosome_size::Int64,
                brkga_params::BrkgaParams, evolutionary_mechanism_on::Bool = true,
    )::BrkgaData

Build a [`BrkgaData`](@ref) object to be used in the evolutionary and path
relink procedures. Such data structure should not be changed outside the
`BrkgaMpIpr` functions. This version accepts all control arguments, and it is
handy for tuning purposes.

# Arguments
- [`problem_instance::AbstractInstance`](@ref AbstractInstance): an instance
  to the problem to be solved.

- `decode_function!::Function`: the decode funtion used to map chromosomes
  to solutions. It **must have** the following signature:

  ```julia
  decode!(chromosome::Array{Float64, 1},
          problem_instance::AbstractInstance,
          rewrite::Bool)::Float64
  ```

  Note that if `rewrite == false`, `decode!` cannot modify `chromosome`.

- [`opt_sense::Sense`](@ref Sense): the optimization sense (
  maximization or minimization).

- `seed::Int64`: seed for the random number generator.

- `chromosome_size::Int64`: number of genes in each chromosome.

- [`brkga_params::BrkgaParams`](@ref BrkgaParams): BRKGA and IPR parameters
  object loaded from a configuration file or manually created. All the data is
  deep-copied.

- `evolutionary_mechanism_on::Bool = true`: if false, no evolution is performed
  but only chromosome decoding. On each generation, all population is replaced
  excluding the best chromosome. Very useful to emulate a multi-start algorithm.

# Throws
- `ArgumentError`: in several cases where the arguments or a combination of them
  are invalid.
"""
function build_brkga(problem_instance::AbstractInstance,
        decode_function!::Function,
        opt_sense::Sense,
        seed::Int64,
        chromosome_size::Int64,
        brkga_params::BrkgaParams,
        evolutionary_mechanism_on::Bool = true,
        )::BrkgaData

    bp = brkga_params

    if evolutionary_mechanism_on
        elite_size = floor(Int64, bp.elite_percentage * bp.population_size)
        num_mutants = floor(Int64, bp.mutants_percentage * bp.population_size)
    else
        elite_size = 1
        num_mutants = bp.population_size - 1
    end

    # Check for errors.
    if chromosome_size < 2
        throw(ArgumentError("chromosome size must be larger than one, " * 
                            "current: $chromosome_size"))
    elseif bp.population_size < 1
        throw(ArgumentError("population size must be larger than zero, " * 
                            "current: $(bp.population_size)"))
    elseif elite_size < 1
        throw(ArgumentError("elite-set size less then one, " * 
                            "current: $elite_size"))
    elseif num_mutants < 0
        throw(ArgumentError("mutant-set size less then zero, " * 
                            "current: $num_mutants"))
    elseif elite_size > bp.population_size
        throw(ArgumentError("elite-set size ($elite_size) greater than " *
                            "population size ($(bp.population_size)"))
    elseif num_mutants >= bp.population_size
        throw(ArgumentError("mutant-set size ($num_mutants) greater than " *
                            "population size ($(bp.population_size)"))
    elseif elite_size + num_mutants > bp.population_size
        throw(ArgumentError("elite ($elite_size) + mutant sets ($num_mutants) " *
                            "greater than population size ($(bp.population_size))"))
    elseif bp.num_elite_parents < 1
        throw(ArgumentError("num_elite_parents must be at least 1, " *
                            "current: $(bp.num_elite_parents)"))
    elseif bp.total_parents < 2
        throw(ArgumentError("total_parents must be at least 2, " *
                            "current: $(bp.total_parents)"))
    elseif bp.num_elite_parents >= bp.total_parents
        throw(ArgumentError("num_elite_parents ($(bp.num_elite_parents) is " *
                            "greater than or equal to " * 
                            "total_parents ($(bp.total_parents)"))
    elseif bp.num_elite_parents > elite_size
        throw(ArgumentError("num_elite_parents ($(bp.num_elite_parents) is " *
                             "greater than elite set ($elite_size)"))
    elseif bp.num_independent_populations < 1
        throw(ArgumentError("number of parallel populations must be larger " *
                            "than zero, current: $(bp.num_independent_populations)"))
    elseif bp.alpha_block_size <= 0.0
        throw(ArgumentError("alpha_block_size must be larger than zero, " *
                            "current: $(bp.alpha_block_size)"))
    elseif bp.pr_percentage <= 0.0 || bp.pr_percentage > 1.0
        throw(ArgumentError("percentage / path size must be in (0, 1], " *
                            "current $(bp.pr_percentage)"))
    end

    brkga_data = BrkgaData(
        opt_sense,
        chromosome_size,
        deepcopy(brkga_params),
        elite_size,
        num_mutants,
        evolutionary_mechanism_on,
        problem_instance,
        decode_function!,
        Random.MersenneTwister(seed),
        Array{Population, 1}(), # previous pop
        Array{Population, 1}(), # current pop
        empty_function,         # bias_function
        0.0,                    # total_bias_weight
        Array{Int64, 1}(undef, bp.population_size),               # shuffled_inds
        Array{Tuple{Float64, Int64}, 1}(undef, bp.total_parents), # parents_ordered
        false,  # initialized
        false   # reset_phase
    )

    if bp.bias_type == LOGINVERSE
        set_bias_custom_function!(brkga_data, r -> 1.0 / log1p(r))

    elseif bp.bias_type == LINEAR
        set_bias_custom_function!(brkga_data, r -> 1.0 / r)

    elseif bp.bias_type == QUADRATIC
        set_bias_custom_function!(brkga_data, r -> r ^ -2.0)

    elseif bp.bias_type == CUBIC
        set_bias_custom_function!(brkga_data, r -> r ^ -3.0)

    elseif bp.bias_type == EXPONENTIAL
        set_bias_custom_function!(brkga_data, r -> exp(-r))

    elseif bp.bias_type == CONSTANT
        set_bias_custom_function!(brkga_data, (::Int64) -> 1.0 / bp.total_parents)
    end

    # Warm up the random number generator.
    rand(brkga_data.rng, 1000)

    return brkga_data
end

################################################################################

"""
    build_brkga(problem_instance, decode_function!, opt_sense, seed,
        chromosome_size, configuration_file,
        evolutionary_mechanism_on)::Tuple{BrkgaData, ExternalControlParams}

Build a [`BrkgaData`](@ref) object to be used in the evolutionary and path
relink procedures, and a [`ExternalControlParams`](@ref) that holds
additional control parameters. Note that [`BrkgaData`](@ref) should not be
changed outside the `BrkgaMpIpr` functions. This version reads most of the
parameters from a configuration file.

# Arguments
- [`problem_instance::AbstractInstance`](@ref AbstractInstance): an instance
  to the problem to be solved.

- `decode_function!::Function`: the decode funtion used to map chromosomes
  to solutions. It **must have** the following signature:

  ```julia
  decode!(chromosome::Array{Float64, 1},
          problem_instance::AbstractInstance,
          rewrite::Bool = true)::Float64
  ```

  Note that if `rewrite == false`, `decode!` cannot modify `chromosome`.

- [`opt_sense::Sense`](@ref Sense): the optimization sense (
  maximization or minimization).

- `seed::Int64`: seed for the random number generator.

- `chromosome_size::Int64`: number of genes in each chromosome..

- `configuration_file::String`:  text file with the BRKGA parameters. An
  example can be found at <a href="example.conf">example.conf</a>. Note that
  the content after "#" is considered comments and it is ignored.

- `evolutionary_mechanism_on::Bool = true`: if false, no evolution is performed
  but only chromosome decoding. On each generation, all population is replaced
  excluding the best chromosome. Very useful to emulate a multi-start algorithm.

# Throws
- `LoadError`: in cases of the file is an invalid configuration file,
  parameters are missing, or parameters are ill-formatted.
- `SystemError`: in case the configuration files cannot be openned.
"""
function build_brkga(
        problem_instance::AbstractInstance,
        decode_function!::Function,
        opt_sense::Sense,
        seed::Int64,
        chromosome_size::Int64,
        configuration_file::String,
        evolutionary_mechanism_on::Bool = true
    )::Tuple{BrkgaData, ExternalControlParams}

    brkga_params, external_params = load_configuration(configuration_file)

    brkga_data = build_brkga(problem_instance, decode_function!, opt_sense,
                             seed, chromosome_size, brkga_params,
                             evolutionary_mechanism_on)

    return (brkga_data, external_params)
end

################################################################################

"""
    initialize!(brkga_data::BrkgaData)

Initialize the populations and others data structures of the BRKGA. If an
initial population is supplied, this method completes the remaining individuals,
if they do not exist.

!!! warning
    THIS METHOD MUST BE CALLED BEFORE ANY OPTIMIZATION METHODS.

This method also performs the initial decoding of the chromosomes. Therefore,
depending on the decoder implementation, this can take a while, and the user may
want to time such procedure in his/her experiments.

!!! note
    As it is in [`evolve!`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using a single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

# Throws
- `ErrorException`: if `bias_function` is not defined previously.
"""
function initialize!(brkga_data::BrkgaData)
    if brkga_data.bias_function == empty_function
        error("the bias function is not defined. Call set_bias_custom_function()!.")
    end

    bd = brkga_data  # Just a short alias.
    params = bd.params

    # If we have warmstaters, complete the population if necessary.
    # Note that it is done only in the true initialization.
    pop_start = 1
    if length(bd.current) > 0 && !bd.reset_phase
        population = bd.current[1]
        for i in (length(population.chromosomes) + 1):params.population_size
            push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
        end
        population.fitness =
            Array{Tuple{Float64, Int64}, 1}(undef, params.population_size)
        pop_start = 2

    elseif length(bd.current) == 0
        bd.current = Array{Population, 1}(undef, params.num_independent_populations)
        bd.previous = Array{Population, 1}(undef, params.num_independent_populations)
    end

    # Build the remaining populations and associated data structures.
    for i in pop_start:params.num_independent_populations
        # If no reset, allocate memory.
        if !bd.reset_phase
            population = Population()
            for j in 1:params.population_size
                push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
            end
            population.fitness =
                Array{Tuple{Float64, Int64}, 1}(undef, params.population_size)
            bd.current[i] = population

        else
            for chr in bd.current[i].chromosomes
                chr .= rand(bd.rng, length(chr))
            end
        end
    end

    # Perform initial decoding. It may take a while.
    for population in bd.current
        Threads.@threads for i in eachindex(population.chromosomes)
            value = bd.decode!(population.chromosomes[i], bd.problem_instance,
                               true)
            population.fitness[i] = (value, i)
        end
        sort!(population.fitness, rev = (bd.opt_sense == MAXIMIZE))
    end

    # Copy the data to previous populations.
    # **NOTE:** (ceandrade) During reset phase, copying item by item maybe
    # faster than deepcoping (which allocates new memory).
    bd.previous = deepcopy(bd.current)

    bd.initialized = true;
    bd.reset_phase = false;
    nothing
end

################################################################################

"""
    set_bias_custom_function!(brkga_data::BrkgaData, bias_function::Function)

Set a new bias function to be used to rank the chromosomes during the mating.
**It must be a positive non-increasing function** returning a Float64, i.e.,
`f(::Int64)::Float64` such that `f(i) ≥ 0` and `f(i) ≥ f(i+1)` for
`i ∈ [1..total_parents]`. For instance, the following sets an inverse quadratic
function:

```julia
set_bias_custom_function!(brkga_data, x -> 1.0 / (x * x))
```

# Throws
- `ArgumentError`: in case the function is not a non-increasing positive
  function.
"""
function set_bias_custom_function!(brkga_data::BrkgaData,
                                   bias_function::Function)

    bias_values = map(bias_function, 1:brkga_data.params.total_parents)

    if any(x -> x < 0.0, bias_values)
        throw(ArgumentError("bias_function must be positive non-increasing"))
    end

    # NOTE (ceandrade) issorted(bias_values, rev=true) and variants do not
    # work for constant functions.
    for i in 2:brkga_data.params.total_parents
        if bias_values[i - 1] < bias_values[i]
           throw(ArgumentError("bias_function is not a non-increasing function"))
        end
    end

    if brkga_data.bias_function != empty_function
        brkga_data.params.bias_type = CUSTOM
    end

    brkga_data.bias_function = bias_function
    brkga_data.total_bias_weight = sum(bias_values)
    nothing
end

################################################################################

"""
    set_initial_population!(brkga_data::BrkgaData,
                            chromosomes::Array{Array{Float64, 1}, 1})

Set initial individuals into the poulation to work as warm-starters. Such
individuals can be obtained from solutions of external procedures such as
fast heuristics, other methaheuristics, or even relaxations from a
mixed integer programming model that models the problem.

All given solutions are assigned to one population only. Therefore, the
maximum number of solutions is the size of the populations.

# Throws
- `ArgumentError`: if the number of given chromosomes is larger than the
  population size; if the sizes of the given chromosomes do not match
  with the required chromosome size.
"""
function set_initial_population!(brkga_data::BrkgaData,
                                 chromosomes::Array{Array{Float64, 1}, 1})

    bd = brkga_data
    if length(chromosomes) > bd.params.population_size
        throw(ArgumentError(
            "number of given chromosomes " *
            "($(length(chromosomes))) is larger than the population size " *
            "($(bd.params.population_size))"
        ))
    end

    # Clean up the current population.
    current = Array{Population, 1}(undef, bd.params.num_independent_populations)
    current[1] = Population()

    for (i, chr) in enumerate(chromosomes)
        if length(chr) != bd.chromosome_size
            msg = "error on setting initial population: chromosome $i does " *
                  "not have the required dimension (actual size: " *
                  "$(length(chr)), required size: $(bd.chromosome_size))"
            throw(ArgumentError(msg))
        end
        push!(current[1].chromosomes, copy(chr))
    end
    bd.current = current
    nothing
end
