################################################################################
# support.jl: initialization, reset, and individual migration routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
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
