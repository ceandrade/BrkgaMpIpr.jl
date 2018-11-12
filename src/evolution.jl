################################################################################
# evolution.jl: main evolutionary routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Apr 19, 2018 by ceandrade
# Last update: Nov 09, 2018 by ceandrade
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
    evolve!(brkga_data::BrkgaData, num_generations::Int64 = 1)

Evolve the current populations following the guidelines of Multi-parent BRKGAs
for `num_generations` generations

# Throws
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `num_generations < 1`.
"""
function evolve!(brkga_data::BrkgaData, num_generations::Int64 = 1)
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if num_generations < 1
        msg = "Number of generations to be evolved must be larger " *
              "than zero: $num_generations"
        throw(ArgumentError(msg))
    end

    for i in 1:num_generations
        for population_index in 1:brkga_data.num_independent_populations
            evolve_population!(brkga_data, population_index)
        end
    end
end

################################################################################

"""
    evolve_population!(brkga_data::BrkgaData, population_index::Int64)

Evolve the population `population_index` to the next.

# Throws
- `ErrorException`: if `initialize!()` was not called before.
- `ArgumentError`: when `population_index < 1` or `population_index > num_independent_populations`.
"""
function evolve_population!(brkga_data::BrkgaData, population_index::Int64)
    bd = brkga_data

    if !bd.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end

    if population_index < 1 || population_index > bd.num_independent_populations
        msg = "Population index must be in [1, " *
              "$(bd.num_independent_populations)]: $population_index"
        throw(ArgumentError(msg))
    end

    # Make names shorter.
    curr = bd.current[population_index]
    next = bd.previous[population_index]

    # First, we copy the elite chromosomes to the next generation.
    for chr in 1:bd.elite_size
        # **NOTE:** the following comment assigments are too slow; it is better
        # do it manually:
        # next.chromosomes[chr] = copy(curr.chromosomes[curr.fitness[chr][2]])
        # next.chromosomes[chr][:] = curr.chromosomes[curr.fitness[chr][2]][:]

        tmp1 = next.chromosomes[chr]
        tmp2 = curr.chromosomes[curr.fitness[chr][2]]
        @inbounds for i in 1:bd.chromosome_size
            tmp1[i] = tmp2[i]
        end

        next.fitness[chr] = (curr.fitness[chr][1], chr);
    end

    # Second, we mate 'pop_size - elite_size - num_mutants' pairs.
    for chr in (bd.elite_size + 1):(bd.population_size - bd.num_mutants)
        # First, we shuffled the elite set and non-elite set indices,
        # then we take the elite and non-elite parents. Note that we cannot
        # shuffled both sets together, otherwise we would mix elite
        # and non-elite individuals.

        bd.shuffled_individuals[1:bd.elite_size] =
            Random.shuffle(bd.rng, 1:bd.elite_size)
        bd.shuffled_individuals[(bd.elite_size + 1):bd.population_size] =
            Random.shuffle(bd.rng, (bd.elite_size + 1):bd.population_size)

        # Take the elite parents.
        @inbounds for i in 1:bd.num_elite_parents
            bd.parents_ordered[i] = curr.fitness[bd.shuffled_individuals[i]]
        end

        # Take the non-elite parents.
        @inbounds for i in 1:(bd.total_parents - bd.num_elite_parents)
            bd.parents_ordered[i + bd.num_elite_parents] =
                curr.fitness[bd.shuffled_individuals[i + bd.elite_size]]
        end

        sort!(bd.parents_ordered, rev = (bd.opt_sense == MAXIMIZE))

        # Performs the mate.
        @inbounds for allele in 1:bd.chromosome_size
            parent = 0
            cumulative_probability = 0.0
            toss = rand(bd.rng)

            while(cumulative_probability < toss)
                parent += 1
                cumulative_probability += bd.bias_function(parent) /
                                          bd.total_bias_weight;
            end

            next.chromosomes[chr][allele] =
                curr.chromosomes[bd.parents_ordered[parent][2]][allele]
        end
    end

    # To finish, we fill up the remaining spots with mutants.
    for chr in (bd.population_size - bd.num_mutants + 1):bd.population_size
        next.chromosomes[chr] .= rand(bd.rng, length(bd.chromosome_size))
    end

    # Perform the decoding on the offpring and mutants.
    Threads.@threads for i in (bd.elite_size + 1):bd.population_size
        value = bd.decode!(next.chromosomes[i], bd.problem_instance)
        next.fitness[i] = (value, i)
    end

    sort!(next.fitness, rev = (bd.opt_sense == MAXIMIZE))

    bd.previous[population_index], bd.current[population_index] =
        bd.current[population_index], bd.previous[population_index]
    nothing
end
