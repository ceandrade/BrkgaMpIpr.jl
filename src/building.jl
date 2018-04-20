################################################################################
# building.jl: Building and parameter setup routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Apr 20, 2018 by ceandrade
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

################################################################################

"""
    build_brkga(problem_instance, decode_function!, opt_sense, seed,
        chromosome_size, population_size, elite_percentage, mutants_percentage,
        evolutionary_mechanism_on, num_elite_parents, total_parents, bias,
        num_independent_populations)::BrkgaData

Build a `BrkgaData` object to be used in the evolutionary and path relink
procedures. Such data structure should not be changed outside the `BrkgaMpIpr`
functions. This version accepts all control arguments, and it is handy for
tuning purposes.

# Arguments
- `problem_instance::AbstractInstance`: an instance to the problem to be
  solved.

- `decode_function!::Function`: the decode funtion used to map chromosomes
  to solutions. It **must have** the following signature:

            decode!(chromosome::Array{Float64, 1},
                    problem_instance::AbstractInstance,
                    rewrite::Bool = true)::Float64

  Note that if `rewrite == false`, `decode!` cannot modify `chromosome`.

- `opt_sense::Sense`: the optimization sense (maximization or minimization).

- `seed::Int64`: seed for the random number generator.

- `chromosome_size::Int64`: number of genes in each chromosome.

- `population_size::Int64`: number of individuals in each population.

- `elite_percentage::Float64`: percentage of elite individuals into each
  population.

- `mutants_percentage::Float64`: percentage of mutants introduced at each
  generation into the population.

- `evolutionary_mechanism_on::Bool = true`: if false, no evolution is performed
  but only chromosome decoding. Very useful to emulate a multi-start algorithm.

- `num_elite_parents::Int64 = 1`: number of parents from elite set used on the
  multi-parent crossover. If `num_elite_parents == 1` and `total_parents == 2`,
  regular BRKGA crossover is performed (one elite and one non-elite sets).

- `total_parents::Int64 = 2`: number of parents to be used on the multi-parent
  crossover. If `total_parents == 2`, regular BRKGA crossover is performed
  (one elite and one non-elite sets).

- `bias::BiasFunction = LOGINVERSE`: bias function used to choose the parents
  during multi-parent crossover. This function substitutes the `rhoe` parameter
  from the original BRKGA.

- `num_independent_populations::Int64 = 1`: number of independent populations
  (island model).
"""
function build_brkga(problem_instance::AbstractInstance,
        decode_function!::Function,
        opt_sense::Sense,
        seed::Int64,
        chromosome_size::Int64,
        population_size::Int64,
        elite_percentage::Float64,
        mutants_percentage::Float64,
        evolutionary_mechanism_on::Bool = true,
        num_elite_parents::Int64 = 1,
        total_parents::Int64 = 2,
        bias::BiasFunction = LOGINVERSE,
        num_independent_populations::Int64 = 1)::BrkgaData

    elite_size = evolutionary_mechanism_on ?
            floor(Int64, elite_percentage * population_size) : 1

    num_mutants = evolutionary_mechanism_on?
            floor(Int64, mutants_percentage * population_size) :
            population_size - 1

    # Check for errors.
    if chromosome_size < 1
        throw(ArgumentError("chromosome size must be larger than zero: $chromosome_size"))
    elseif population_size < 1
        throw(ArgumentError("population size must be larger than zero: $population_size"))
    elseif elite_size < 1
        throw(ArgumentError("elite-set size less then one: $elite_size"))
    elseif num_mutants < 0
        throw(ArgumentError("mutant-set size less then zero: $num_mutants"))
    elseif elite_size > population_size
        throw(ArgumentError("elite-set size ($elite_size) greater than population size ($population_size)"))
    elseif num_mutants >= population_size
        throw(ArgumentError("mutant-set size ($num_mutants) greater than population size ($population_size)"))
    elseif elite_size + num_mutants > population_size
        throw(ArgumentError("elite ($elite_size) + mutant sets ($num_mutants) greater than population size ($population_size)"))
    elseif num_elite_parents < 1
        throw(ArgumentError("num_elite_parents must be at least 1: $num_elite_parents"))
    elseif total_parents < 2
        throw(ArgumentError("total_parents must be at least 2: $total_parents"))
    elseif num_elite_parents >= total_parents
        throw(ArgumentError("num_elite_parents ($num_elite_parents) is greater than or equal to total_parents ($total_parents)"))
    elseif num_elite_parents > elite_size
        throw(ArgumentError("num_elite_parents ($num_elite_parents) is greater than elite set ($elite_size)"))
    elseif num_independent_populations < 1
        throw(ArgumentError("number of parallel populations must be larger than zero: $num_independent_populations"))

    # TODO (ceandrade): check path relink params here
    end

    brkga_data = BrkgaData(
        opt_sense,
        chromosome_size,
        population_size,
        elite_size,
        num_mutants,
        num_elite_parents,
        total_parents,
        bias,
        num_independent_populations,
        evolutionary_mechanism_on,
        # TODO (ceandrade): list the IPR parameters here.
        problem_instance,
        decode_function!,
        MersenneTwister(seed),
        Array{Population, 1}(),                            # previous pop
        Array{Population, 1}(),                            # current pop
        empty_function,                                    # bias_function
        0.0,                                               # total_bias_weight
        Array{Int64, 1}(population_size),                  # shuffled_inds
        Array{Tuple{Float64, Int64}, 1}(total_parents),    # parents_ordered
        false,                                             # initialized
        false                                              # reset_phase
    )

    if bias == LOGINVERSE
        set_bias_custom_function!(brkga_data, r -> 1.0 / log1p(r))

    elseif bias == LINEAR
        set_bias_custom_function!(brkga_data, r -> 1.0 / r)

    elseif bias == QUADRATIC
        set_bias_custom_function!(brkga_data, r -> r ^ -2.0)

    elseif bias == CUBIC
        set_bias_custom_function!(brkga_data, r -> r ^ -3.0)

    elseif bias == EXPONENTIAL
        set_bias_custom_function!(brkga_data, r -> exp(-r))

    elseif bias == CONSTANT
        set_bias_custom_function!(brkga_data, (::Int64) -> 1.0 / total_parents)
    end

    # Warm up the random number generator.
    rand(brkga_data.rng, 1000)

    return brkga_data
end

"""
    build_brkga(problem_instance, decode_function!, opt_sense, seed,
        chromosome_size, configuration_file,
        evolutionary_mechanism_on)::Tuple{BrkgaData, ExternalControlParams}

Build a `BrkgaData` object to be used in the evolutionary and path relink
procedures, and a `ExternalControlParams` that holds additional control
parameters. Note that `BrkgaData` should not be changed outside the `BrkgaMpIpr`
functions. This version reads most of the parameters from a configuration file.

# Arguments
- `problem_instance::AbstractInstance`: an instance to the problem to be
  solved.

- `decode_function!::Function`: the decode funtion used to map chromosomes
  to solutions. It **must have** the following signature:

            decode!(chromosome::Array{Float64, 1},
                    problem_instance::AbstractInstance,
                    rewrite::Bool = true)::Float64

  Note that if `rewrite == false`, `decode!` cannot modify `chromosome`.

- `opt_sense::Sense`: the optimization sense (maximization or minimization).

- `seed::Int64`: seed for the random number generator.

- `chromosome_size::Int64`: number of genes in each chromosome.

- `configuration_file::String`:  text file with the BRKGA parameters. An
  example can be found at <a href="example.conf">example.conf</a>. Note that
  the content after "#" is considered comments and it is ignored.

 - `evolutionary_mechanism_on::Bool = true`: if false, no evolution is performed
   but only chromosome decoding. Very useful to emulate a multi-start algorithm.

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

    # TODO (ceandrade) add the path relink parameters.

    const param_names_types = [
        ("POPULATION_SIZE", Int64),
        ("ELITE_PERCENTAGE", Float64),
        ("MUTANTS_PERCENTAGE", Float64),
        ("ELITE_PARENTS", Int64),
        ("TOTAL_PARENTS", Int64),
        ("BIAS_FUNCTION", BiasFunction),
        ("INDEPENDENT_POPULATIONS", Int64),
        # ("PR_MINIMUM_DISTANCE", Float64),
        # ("PR_TYPE", ),
        # ("PR_SELECTION", ),
        # ("ALPHA_BLOCK_SIZE", Float64),
        # ("PR_PERCENTAGE", Float64),
        ("EXCHANGE_INTERVAL", Int64),
        ("NUM_EXCHANGE_INDIVUDUALS", Int64),
        ("RESET_INTERVAL", Int64),
    ]

    const param_index = Dict(
        [v[1] => i for (i, v) in enumerate(param_names_types)])

    param_values = Array{Any}(length(param_names_types))
    param_given = falses(length(param_names_types))

    lines = Array{String,1}()
    open(configuration_file) do file
        lines = readlines(file)
    end
    if length(lines) == 0
        throw(LoadError(configuration_file, 0, "cannot read '$configuration_file'"))
    end

    for (line_number, line) in enumerate(lines)
        line = uppercase(strip(line))
        if length(line) == 0 || line[1] == '#'
            continue
        end

        param_name = ""
        value = 0
        try
            param_name, value = split(line)
            idx = param_index[param_name]
            param_values[idx] = parse(param_names_types[idx][2], String(value))
            param_given[idx] = true
        catch err
            if isa(err, BoundsError)
                throw(LoadError(configuration_file, line_number,
                        "error line $line_number of '$configuration_file': missing parameter or value"))

            elseif isa(err, KeyError)
                throw(LoadError(configuration_file, line_number, "parameter '$param_name' unknown"))

            elseif isa(err, ArgumentError)
                throw(LoadError(configuration_file, line_number, "invalid value for '$param_name': $value"))
            end
        end
    end

    missing_params = ""
    for (idx, value) in enumerate(param_given)
        if !value
            missing_params *= "'" * param_names_types[idx][1] * "',"
        end
    end
    if length(missing_params) > 0
        throw(LoadError(configuration_file, 0, "missing parameters: $missing_params"))
    end

    external_params = ExternalControlParams(
        param_values[param_index["EXCHANGE_INTERVAL"]],
        param_values[param_index["NUM_EXCHANGE_INDIVUDUALS"]],
        param_values[param_index["RESET_INTERVAL"]]
    )

    brkga_data = build_brkga(problem_instance, decode_function!, opt_sense,
                    seed, chromosome_size,
                    param_values[param_index["POPULATION_SIZE"]],
                    param_values[param_index["ELITE_PERCENTAGE"]],
                    param_values[param_index["MUTANTS_PERCENTAGE"]],
                    evolutionary_mechanism_on,
                    param_values[param_index["ELITE_PARENTS"]],
                    param_values[param_index["TOTAL_PARENTS"]],
                    param_values[param_index["BIAS_FUNCTION"]],
                    param_values[param_index["INDEPENDENT_POPULATIONS"]])

    return (brkga_data, external_params)
end

################################################################################

"""
    initialize!(brkga_data::BrkgaData)

Initialize the populations and others data structures of the BRKGA. If an
initial population is supplied, this method completes the remaining individuals,
if they do not exist.

**THIS METHOD MUST BE CALLED BEFORE ANY OPTIMIZATIOM METHODS.**

This method also performs the initial decoding of the chromosomes. Therefore,
depending on the decoder implementation, this can take a while, and the user may
want to time such procedure in his/her experiments.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).

# Throws
- `ErrorException`: if `bias_function` is not defined previously.
"""
function initialize!(brkga_data::BrkgaData)
    if brkga_data.bias_function == empty_function
        error("bias function is not defined. Call set_bias_custom_function()!.")
    end

    bd = brkga_data  # Just a short alias.

    # If we have warmstaters, complete the population if necessary.
    # Note that it is done only in the true initialization.
    pop_start = 1
    if length(bd.current) > 0 && !bd.reset_phase
        population = bd.current[1]
        for i = (length(population.chromosomes) + 1):bd.population_size
            push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
        end
        population.fitness = Array{Tuple{Float64, Int64}, 1}(bd.population_size)
        pop_start = 2

    elseif length(bd.current) == 0
        bd.current = Array{Population, 1}(bd.num_independent_populations)
        bd.previous = Array{Population, 1}(bd.num_independent_populations)
    end

    # Build the remaining populations and associated data structures.
    for i = pop_start:bd.num_independent_populations
        # If no reset, allocate memory.
        if !bd.reset_phase
            population = Population()
            for j = 1:bd.population_size
                push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
            end
            population.fitness =
                Array{Tuple{Float64, Int64}, 1}(bd.population_size)
            brkga_data.current[i] = population

        else
            for chr in brkga_data.current[i].chromosomes
                chr .= rand(bd.rng, length(chr))
            end
        end
    end

    # Perform initial decoding. It may take a while.
    for population in bd.current
        Threads.@threads for i in eachindex(population.chromosomes)
            value = bd.decode!(population.chromosomes[i], bd.problem_instance)
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
**It must be a positive non-decreasing function** returning a Float64, i.e.,
`f(::Int64)::Float64` such that `f(i) ≥ 0` and `f(i) ≥ f(i+1)` for
`i ∈ [1..total_parents]`.

# Throws
- `ArgumentError`: in case the function is not a non-decreasing positive
  function.
"""
function set_bias_custom_function!(brkga_data::BrkgaData,
                                   bias_function::Function)

    bias_values = map(bias_function, 1:brkga_data.total_parents)
    if any(x -> x < 0.0, bias_values)
        throw(ArgumentError("bias_function must be positive non-decreasing"))
    end

    # NOTE (ceandrade) issorted(bias_values, rev=true) and variants do not
    # work for constant functions.
    for i in 2:brkga_data.total_parents
        if bias_values[i - 1] < bias_values[i]
           throw(ArgumentError("bias_function is not a non-decreasing function"))
        end
    end

    if brkga_data.bias_function != empty_function
        brkga_data.bias = CUSTOM
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

# Throws
- `ArgumentError`: if the number of given chromosomes is larger than the
  population size; if the sizes of the given chromosomes do not match
  with the required chromosome size.
"""
function set_initial_population!(brkga_data::BrkgaData,
                                 chromosomes::Array{Array{Float64, 1}, 1})

    if length(chromosomes) > brkga_data.population_size
        throw(ArgumentError("number of given chromosomes " *
            "($(length(chromosomes))) is larger than the population size " *
            "($(brkga_data.population_size))"))
    end

    # Clean up the current population.
    current = Array{Population, 1}(brkga_data.num_independent_populations)
    current[1] = Population()

    for (i, chr) in enumerate(chromosomes)
        if length(chr) != brkga_data.chromosome_size
            msg = "error on setting initial population: chromosome $i does " *
                "not have the required dimension (actual size: " *
                "$(length(chr)), required size: $(brkga_data.chromosome_size))"
            throw(ArgumentError(msg))
        end
        push!(current[1].chromosomes, copy(chr))
    end
    brkga_data.current = current
    nothing
end
