################################################################################
# support.jl: initialization, reset, and individual migration routines.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
# Last update: Feb 08, 2019 by ceandrade
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
    reset!(brkga_data::BrkgaData)

Reset all populations with brand new keys. All warm-start solutions provided
by [`set_initial_population!`](@ref) are discarded. You may use
[`inject_chromosome!`](@ref) to insert those solutions again.

!!! warning
    As it is in [`evolve!()`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/)

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
"""
function reset!(brkga_data::BrkgaData)
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before reset!()")
    end
    brkga_data.reset_phase = true
    initialize!(brkga_data)
    nothing
end

################################################################################

"""
    exchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)

Exchange elite-solutions between the populations. Given a population, the
`num_immigrants` best solutions are copied to the neighbor populations,
replacing their worth solutions. If there is only one population, nothing is
done.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
- `ArgumentError`: when `num_immigrants < 1`.
- `ArgumentError`: `num_immigrants ≥ ⌈population_size/num_independent_populations⌉ - 1`.
"""
function exchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)
    bd = brkga_data  # Just a short alias.

    if !bd.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before exchange_elite!()")
    end

    if bd.params.num_independent_populations == 1
        return nothing
    end

    immigrants_threshold = cld(bd.params.population_size,
                               bd.params.num_independent_populations - 1)

    if num_immigrants < 1 || num_immigrants >= immigrants_threshold
        msg = "number of immigrants ($num_immigrants) less than one, or " *
              "larger than or equal to population size / " *
              "num_independent_populations ($immigrants_threshold)"
        throw(ArgumentError(msg))
    end

    # Population i receives num_immigrants best individuals from j
    # overwriting worst i individuals.
    @inbounds for i in 1:bd.params.num_independent_populations
        to_pop = bd.current[i]
        dest = bd.params.population_size
        @inbounds for j in 1:bd.params.num_independent_populations
            if j == i
                continue
            end

            @inbounds for m in 1:num_immigrants
                from_pop = bd.current[j]
                (from_value, from_idx) = from_pop.fitness[m]
                to_idx = to_pop.fitness[dest][2]

                # Copy keys.
                to_pop.chromosomes[to_idx] .= from_pop.chromosomes[from_idx]
                # Keep the same index but with new value.
                to_pop.fitness[dest] = (from_value, to_idx)
                dest -= 1
            end
        end
    end

    # Resort.
    Threads.@threads for i in 1:bd.params.num_independent_populations
        sort!(bd.current[i].fitness, rev = bd.opt_sense == MAXIMIZE)
    end
    nothing
end

################################################################################

"""
    get_best_fitness!(brkga_data::BrkgaData)::Float64

Return the fitness/value of the best individual found so far among all
populations.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
"""
function get_best_fitness(brkga_data::BrkgaData)::Float64
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before get_best_fitness!()")
    end

    best = brkga_data.current[1].fitness[1][1]
    for i in 2:brkga_data.params.num_independent_populations
        if (brkga_data.current[i].fitness[1][1] < best) ==
           (brkga_data.opt_sense == MINIMIZE)
            best = brkga_data.current[i].fitness[1][1]
        end
    end
    return best
end

################################################################################

"""
    get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}

Return a copy of the best individual found so far among all populations.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
"""
function get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before get_best_chromosome!()")
    end

    best_population = 1
    best_individual = brkga_data.current[1].fitness[1][2]
    for i in 2:brkga_data.params.num_independent_populations
        fitness, individual = brkga_data.current[i].fitness[1]
        if (fitness < brkga_data.current[best_population].fitness[1][1]) ==
           (brkga_data.opt_sense == MINIMIZE)
            best_population = i
            best_individual = individual
        end
    end

    return copy(brkga_data.current[best_population].
                chromosomes[best_individual])
end

################################################################################

"""
    get_chromosome(brkga_data::BrkgaData, population_index::Int64,
                   position::Int64)::Array{Float64, 1}

Return a copy of the chromosome `position` in the population `population_index`.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
- `ArgumentError`: when `population_index < 1` or
  `population_index > num_independent_populations`.
- `ArgumentError`: when `position < 1` or `position > population_size`.
"""
function get_chromosome(brkga_data::BrkgaData, population_index::Int64,
                        position::Int64)::Array{Float64, 1}
    bd = brkga_data
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before get_chromosome!()")
    end

    if population_index < 1 ||
       population_index > bd.params.num_independent_populations
        msg = "Population must be in [1, " *
              "$(bd.params.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    if position < 1 || position > bd.params.population_size
        msg = "Chromosome position must be in [1, " *
              "$(bd.params.population_size)]: $position"
        throw(ArgumentError(msg))
    end

    pop = bd.current[population_index]
    return copy(pop.chromosomes[pop.fitness[position][2]])
end

################################################################################

"""
    get_current_population(brkga_data::BrkgaData,
                           population_index::Int64)::Population

Return a reference for population `population_index`.

!!! note
    This function is implemented for complaince with the C++ API. The user
    can access the population directly using
    `brkga_data.current[population_index]`.

!!! warning
    IT IS NOT ADIVISED TO CHANGE THE POPULATION DIRECTLY, since such changes
    can result in undefined behavior.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
- `ArgumentError`: when `population_index < 1` or
  `population_index > num_independent_populations`.
"""
function get_current_population(brkga_data::BrkgaData,
                                population_index::Int64)::Population
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before get_current_population!()")
    end

    if population_index < 1 ||
       population_index > brkga_data.params.num_independent_populations
        msg = "Population must be in [1, " *
              "$(brkga_data.params.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    return brkga_data.current[population_index]
end

################################################################################

"""
    function inject_chromosome!(brkga_data::BrkgaData,
                                chromosome::Array{Float64, 1},
                                population_index::Int64,
                                position::Int64,
                                fitness::Float64 = Inf)

Inject `chromosome` and its `fitness` into population `population_index`
in the `position` place. If fitness is not provided (`fitness = Inf`), the
decoding is performed over `chromosome`. Once the chromosome is injected,
the population is re-sorted according to the chromosomes' fitness.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
- `ArgumentError`: when `population_index < 1` or
  `population_index > num_independent_populations`.
- `ArgumentError`: when `position < 1` or `position > population_size`.
- `ArgumentError`: when `lenght(chromosome) != chromosome_size`.
"""
function inject_chromosome!(brkga_data::BrkgaData,
                            chromosome::Array{Float64, 1},
                            population_index::Int64,
                            position::Int64,
                            fitness::Float64 = Inf)
    bd = brkga_data
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before inject_chromosome!()")
    end

    if population_index < 1 ||
       population_index > bd.params.num_independent_populations
        msg = "Population must be in [1, " *
              "$(bd.params.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    if position < 1 || position > bd.params.population_size
        msg = "Chromosome position must be in [1, " *
              "$(bd.params.population_size)]: $position"
        throw(ArgumentError(msg))
    end

    if length(chromosome) != bd.chromosome_size
        msg = "Chromosome size must be $(bd.chromosome_size) not " *
              "$(length(chromosome))"
        throw(ArgumentError(msg))
    end

    pop = bd.current[population_index]
    idx = pop.fitness[position][2]
    pop.chromosomes[idx] .= chromosome

    if fitness == Inf
        fitness = bd.decode!(pop.chromosomes[idx], bd.problem_instance)
    end

    pop.fitness[position] = (fitness, idx)
    sort!(pop.fitness, rev = (bd.opt_sense == MAXIMIZE))
    nothing
end

################################################################################

"""
    function shake!(brkga_data::BrkgaData, intensity::Int64,
                    shaking_type::ShakingType, population_index::Int64 = Inf64)

Perform a shaking in the chosen population. The procedure applies changes
(shaking) on elite chromosomes and fully reset the remaining population.

# Arguments
- `intensity::Int64`: the intensity of the shaking (> 0);
- [`shaking_type::ShakingType`](@ref ShakingType): either `CHANGE` or `SWAP`
  moves;
- `population_index::Int64`: the index of the population to be shaken.
  If `population_index > num_independent_populations`, all populations are
  shaken.

# Throws
- `ErrorException`: if [`initialize!`](@ref) has not been called before.
- `ArgumentError`: when `population_index < 1` or `intensity < 1`.
"""
function shake!(brkga_data::BrkgaData, intensity::Int64,
                shaking_type::ShakingType, population_index::Int64 = Inf64)
    bd = brkga_data
    if !bd.initialized
        error("the algorithm hasn't been initialized. Call initialize!() " *
              "before shake!()")
    end

    if intensity < 1
        msg = "intensity must be larger than zero: $intensity"
        throw(ArgumentError(msg))
    end

    if population_index < 1
        msg = "population index must be larger than zero: $population_index"
        throw(ArgumentError(msg))
    end

    if population_index > bd.params.num_independent_populations
        pop_start = 1
        pop_end = bd.params.num_independent_populations
    else
        pop_start = population_index
        pop_end = population_index
    end

    for pop_idx in pop_start:pop_end
        pop = brkga_data.current[pop_idx].chromosomes

        # Shake the elite set.
        for e in 1:bd.elite_size, k in 1:intensity
            i = rand(bd.rng, 1:(bd.chromosome_size - 1))
            if shaking_type == CHANGE
                # Invert value.
                pop[e][i] = 1.0 - pop[e][i]
            else
                # Swap with neighbor.
                pop[e][i], pop[e][i + 1] = pop[e][i + 1], pop[e][i]
            end

            i = rand(bd.rng, 1:bd.chromosome_size)
            if shaking_type == CHANGE
                # Change to random value.
                pop[e][i] = rand(bd.rng)
            else
                # Swap two random positions.
                j = rand(bd.rng, 1:bd.chromosome_size)
                pop[e][i], pop[e][j] = pop[e][j], pop[e][i]
            end
        end

        # Reset the remaining population.
        for chr in (bd.elite_size + 1):bd.params.population_size
            pop[chr] .= rand(bd.rng, bd.chromosome_size)
        end

        # Perform the decoding.
        curr = brkga_data.current[pop_idx]
        Threads.@threads for i in 1:bd.params.population_size
            value = bd.decode!(curr.chromosomes[i], bd.problem_instance)
            curr.fitness[i] = (value, i)
        end

        sort!(curr.fitness, rev = (bd.opt_sense == MAXIMIZE))
    end
    nothing
end
