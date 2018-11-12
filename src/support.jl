################################################################################
# support.jl: initialization, reset, and individual migration routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
# Last update: Nov 08, 2018 by ceandrade
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

Reset all populations with brand new keys. All warm start solutions provided
by `set_initial_population!()` are discarded.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).

# Throws
- `ErrorException`: if `initialize!()` was not called before.
"""
function reset!(brkga_data::BrkgaData)
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
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
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `num_immigrants < 1`.
- `ArgumentError`: `num_immigrants ≥ ⌈population_size/num_independent_populations⌉ - 1`.
"""
function exchange_elite!(brkga_data::BrkgaData, num_immigrants::Int64)
    bd = brkga_data  # Just a short alias.

    if !bd.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if brkga_data.num_independent_populations == 1
        return nothing
    end

    immigrants_threshold = cld(bd.population_size,
                               bd.num_independent_populations - 1)

    if num_immigrants < 1 || num_immigrants >= immigrants_threshold
        msg = "number of immigrants ($num_immigrants) less than one, or " *
              "larger than or equal to population size / " *
              "num_independent_populations ($immigrants_threshold)"
        throw(ArgumentError(msg))
    end

    # Population i receives num_immigrants best individuals from j
    # overwriting worst i individuals.
    for i = 1:bd.num_independent_populations
        to_pop = bd.current[i]
        dest = bd.population_size
        for j = 1:bd.num_independent_populations
            if j == i
                continue
            end

            for m = 1:num_immigrants
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
    Threads.@threads for i = 1:bd.num_independent_populations
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
- `ErrorException`: if `initialize!()` was not called before.
"""
function get_best_fitness(brkga_data::BrkgaData)::Float64
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    best = brkga_data.current[1].fitness[1][1]
    for i in 2:brkga_data.num_independent_populations
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
- `ErrorException`: if `initialize!()` was not called before.
"""
function get_best_chromosome(brkga_data::BrkgaData)::Array{Float64, 1}
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    best_population = 1
    best_individual = 1
    for i in 2:brkga_data.num_independent_populations
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
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `population_index < 1` or `population_index > num_independent_populations`.
- `ArgumentError`: when `position < 1` or `position > population_size`.
"""
function get_chromosome(brkga_data::BrkgaData, population_index::Int64,
                        position::Int64)::Array{Float64, 1}
    bd = brkga_data
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if population_index < 1 || population_index > bd.num_independent_populations
        msg = "Population must be in [1, " *
              "$(bd.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    if position < 1 || position > bd.population_size
        msg = "Chromosome position must be in [1, " *
              "$(bd.population_size)]: $position"
        throw(ArgumentError(msg))
    end

    pop = bd.current[population_index]
    return copy(pop.chromosomes[pop.fitness[position][2]])
end

################################################################################

"""
    get_current_population(brkga_data::BrkgaData,
                           population_index::Int64)::Population

Return a reference for population `population_index`. **NOTE:** this function is
implemented for complaince with the C++ API. The user can access the population
directly using `brkga_data.current[population_index]`. However, **IT IS NOT
ADIVISED TO CHANGE THE POPULATION DIRECTLY,** since such changes can result in
undefined behavior.

# Throws
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `population_index < 1` or `population_index > num_independent_populations`.
"""
function get_current_population(brkga_data::BrkgaData,
                                population_index::Int64)::Population
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if population_index < 1 ||
       population_index > brkga_data.num_independent_populations
        msg = "Population must be in [1, " *
              "$(brkga_data.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    return brkga_data.current[population_index]
end

################################################################################

"""
    function inject_chromosome(brkga_data::BrkgaData,
                               chromosome::Array{Float64, 1},
                               population_index::Int64,
                               position::Int64,
                               fitness::Float64 = Inf)

Inject `chromosome` and its `fitness` into population `population_index`
in the `position` place. If fitness is not provided (`fitness = Inf`), the
decoding is performed over `chromosome`. Once the chromosome is injected,
the population is re-sorted according to the chromosomes' fitness.

# Throws
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `population_index < 1` or `population_index > num_independent_populations`.
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
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if population_index < 1 || population_index > bd.num_independent_populations
        msg = "Population must be in [1, " *
              "$(bd.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    if position < 1 || position > bd.population_size
        msg = "Chromosome position must be in [1, " *
              "$(bd.population_size)]: $position"
        throw(ArgumentError(msg))
    end

    if length(chromosome) != bd.chromosome_size
        msg = "Chromosome size must be $(bd.chromosome_size) not " *
              "$(length(chromosome))"
        throw(ArgumentError(msg))
    end

    pop = bd.current[population_index]
    idx = pop.fitness[position][2]
    unsafe_copyto!(pop.chromosomes[idx], 1, chromosome, 1, bd.chromosome_size)

    if fitness == Inf
        fitness = bd.decode!(pop.chromosomes[idx], bd.problem_instance)
    end

    pop.fitness[position] = (fitness, idx)
    sort!(pop.fitness, rev = (bd.opt_sense == MAXIMIZE))
    nothing
end
