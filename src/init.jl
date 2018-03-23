################################################################################
# init.jl: Building andiInitialization routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Mar 23, 2018 by ceandrade
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
    build_brkga(problem_instance, decode_function!, opt_sense, seed,
         chromosome_size, population_size, elite_percentage, mutants_percentage,
         evolutionary_mechanism_on, num_elite_parents, total_parents, bias,
         num_independent_populations)::BrkgaData

Build a `BrkgaData` object to be used in the evolutionary and path relink
procedures. Such data structure should not be changed outside the `BrkgaMpIpr`
functions. This version accepts all control arguments, and it is handy for
tuning purposes.

# Arguments:
- `problem_instance::AbstractInstance`: an instance to the problem to be
  solved.

- `decode_function!::Function`: the decode funtion used to map chromosomes
  to solutions. It **must have** the following signature:

             decode!(chromosome::Array{Float64, 1},
                     problem_instance::AbstractInstance,
                     rewrite::Bool = true)::Float64

  Note that if `rewrite == true`, `decode!` can change `chromosome`.

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
        Array{Population, 1}(num_independent_populations),
        Array{Population, 1}(num_independent_populations),
        x -> 0.1,
        0.0,
        Array{Int64, 1}(population_size),
        Array{Tuple{Float64,Int64}, 1}(population_size),
        false,
        false,
        false
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
        set_bias_custom_function!(brkga_data, r -> 1.0 / total_parents)
    end

    # Warm up the random number generator.
    rand(brkga_data.rng, 1000)

    return brkga_data
end

# """
#     build_brkga(problem_instance, decode_function!, opt_sense, seed,
#          chromosome_size, configuration_file,
#          evolutionary_mechanism_on)::BrkgaData
#
# Build a `BrkgaData` object to be used in the evolutionary and path relink
# procedures. Such data structure should not be changed outside the `BrkgaMpIpr`
# functions. This version reads most of the parameters from a configuration file.
#
# # Arguments:
# - `problem_instance::AbstractInstance`: an instance to the problem to be
#   solved.
#
# - `decode_function!::Function`: the decode funtion used to map chromosomes
#   to solutions. It **must have** the following signature:
#
#              decode!(chromosome::Array{Float64, 1},
#                      problem_instance::AbstractInstance,
#                      rewrite::Bool = true)::Float64
#
#   Note that if `rewrite == true`, `decode!` can change `chromosome`.
#
# - `opt_sense::Sense`: the optimization sense (maximization or minimization).
#
# - `seed::Int64`: seed for the random number generator.
#
# - `chromosome_size::Int64`: number of genes in each chromosome.
#
# - `configuration_file::String`:  text file with the BRKGA parameters. An
#   example can be found at <a href="example.conf">example.conf</a>. Note that
#   the content after "#" is considered comments and it is ignored.
#
#  - `evolutionary_mechanism_on::Bool = true`: if false, no evolution is performed
#    but only chromosome decoding. Very useful to emulate a multi-start algorithm.
# """
# function build_brkga(problem_instance::AbstractInstance,
#               decode_function!::Function,
#               opt_sense::Sense,
#               seed::Int64,
#               chromosome_size::Int64,
#               configuration_file::String,
#               evolutionary_mechanism_on::Bool = true)::BrkgaData
#
#     brkga_data = BrkgaData()
#     # brkga_data.opt_sense = opt_sense
#     #
#     # brkga_data.rng = MersenneTwister(1234)
#
#     # brkga_data.num_chromosomes = num_chromosomes
#     # brkga_data.chromosome_size = chromosome_size
#     #
#     # brkga_data.chromosomes = Array{Array{Float64}}(0)
#     # for i in 1:num_chromosomes
#     #     push!(brkga_data.chromosomes, rand(brkga_data.rng, chromosome_size))
#     # end
#
#     return brkga_data
# end

################################################################################

"""
    set_bias_custom_function!(brkga_data::BrkgaData, bias_function::Function)

Set a new bias function to be used to rank the chromosomes during the mating.
**It must be a positive non-decreasing function** returning a Float64, i.e.,
`f(::Int64)::Float64` such that `f(i) ≥ 0` and `f(i) ≥ f(i+1)` for
`i ∈ [1..total_parents]`.
"""
function set_bias_custom_function!(brkga_data::BrkgaData,
                                   bias_function::Function)

    bias_values = map(bias_function, 1:brkga_data.total_parents)

    if any(x -> x < 0.0, bias_values)
        throw(ArgumentError("bias_function must be positive non-decreasing"))
    end

    # TODO (ceandrade) issorted(bias_values, rev=true) and variants do not
    # work for constant functions.
    for i in 2:brkga_data.total_parents
         if bias_values[i - 1] < bias_values[i]
             throw(ArgumentError("bias_function is not a non-decreasing function"))
         end
    end

    brkga_data.bias_function = bias_function
    brkga_data.total_bias_weight = sum(bias_values)
end
