################################################################################
# support.jl: initialization, reset, and individual migration routines.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
# Last update: Feb 28, 2019 by ceandrade
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

Return a copy of the chromosome ranked at `position` in the population
`population_index`. Note that the chromosomes are rakend by fitness and
the best chromosome is located in position 1.

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
