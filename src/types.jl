################################################################################
# types.jl: Definitions of internal data structures and external API.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Feb 12, 2019 by ceandrade
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
# Enumerations
################################################################################

"""
    @enum Sense

Tells the algorithm either to `MINIMIZE` or `MAXIMIZE` the objective function.
"""
@enum Sense begin
   MINIMIZE = 0
   MAXIMIZE = 1
end

################################################################################

"""
    @enum BiasFunction

Specifies a bias function when choosing parents to mating. This function
substitutes the ``\\rho`` (rho) parameter from the original BRKGA.
For a given rank ``r``, we have the following functions:

  - `CONSTANT`: 1 / number of parents for mating (all individuals have the
                same probability)

  - `CUBIC`: ``r^{-3}``

  - `EXPONENTIAL`: ``ϵ^{-r}``

  - `LINEAR`: ``1 / r``

  - `LOGINVERSE`: ``1 / \\log(r + 1)``

  - `QUADRATIC`: ``r^{-2}``
"""
@enum BiasFunction begin
    CONSTANT = 0
    CUBIC = 1
    EXPONENTIAL = 2
    LINEAR = 3
    LOGINVERSE = 4
    QUADRATIC = 5
    CUSTOM = 6
end

################################################################################

"""
    @enum PathRelinkingType

Specifies type of path relinking:

  - `DIRECT`: changes each key for the correspondent in the other chromosome.

  - `PERMUTATION`: switches the order of a key for that in the other chromosome.
"""
@enum PathRelinkingType begin
    DIRECT = 0
    PERMUTATION = 1
end

################################################################################

"""
    @enum PathRelinkingSelection

Specifies which individuals used to build the path:

  - `BESTSOLUTION`: selects, in the order, the best solution of each population.

  - `RANDOMELITE`: chooses uniformly random solutions from the elite sets.
"""
@enum PathRelinkingSelection begin
    BESTSOLUTION = 0
    RANDOMELITE = 1
end

################################################################################

"""
    @enum PathRelinkingResult

Specifies the result type/status of path relink procedure:

  - `TOO_HOMOGENEOUS`: the chromosomes among the populations are too homogeneous
                       and the path relink will not generate improveded
                       solutions.

  - `NO_IMPROVEMENT`: path relink was done but no improveded solution was found.

  - `ELITE_IMPROVEMENT`: an improved solution among the elite set was found,
                         but the best solution was not improved.

  - `BEST_IMPROVEMENT`: the best solution was improved.
"""
@enum PathRelinkingResult begin
    TOO_HOMOGENEOUS = 0
    NO_IMPROVEMENT = 1
    ELITE_IMPROVEMENT = 3
    BEST_IMPROVEMENT = 7
end

################################################################################

"""
    @enum ShakingType

Specifies the type of shaking to be performed.

  - `CHANGE`: applies the following perturbations:
    1) Inverts the value of a random chosen, i.e., from `value` to
       `1 - value`;
    2) Assigns a random value to a random key.

  - `SWAP`: applies two swap perturbations:
    1) Swaps the values of a randomly chosen key `i` and its
       neighbor `i + 1`;
    2) Swaps values of two randomly chosen keys.
"""
@enum ShakingType begin
    CHANGE = 0
    SWAP = 1
end

################################################################################
# Data structures
################################################################################

"""
    abstract type AbstractInstance

The required interface for external data to be provided to the decoder. The
problem definitions and data **must be** a subtype of AbstractInstance. For
example,

```julia
mutable struct TSPInstance <: AbstractInstance
    num_cities::Int64
    distances::Array{Float64}
end
```
represents an instance type for the Traveling Salesman problem which defines
the number fo cities and a matrix of distances between them.
"""
abstract type AbstractInstance end

################################################################################

"""
    mutable struct Population (internal BRKGA data struct)

Encapsulates a population of chromosomes. Note that this struct is **NOT**
meant to be used externally of this unit.

Fields
------
$(FIELDS)
"""
mutable struct Population
    """
    Population of chromosomes.
    """
    chromosomes::Array{Array{Float64, 1}, 1}

    """
    Fitness of a each chromosome. Each pair represents the fitness and
    the chromosome index.
    """
    fitness::Array{Tuple{Float64, Int64}, 1}

    """
    Default constructor.
    """
    Population() = new(Array{Array{Float64, 1}, 1}(),
                       Array{Tuple{Float64, Int64}, 1}())

    """
    Copy constructor.
    """
    Population(pop::Population) = new(deepcopy(pop.chromosomes),
                                      deepcopy(pop.fitness))
end

################################################################################

"""
    mutable struct BrkgaParams

Represents the BRKGA and IPR hyper-parameters. You can load these parameters
from a configuration file using [`load_configuration`](@ref) and
[`build_brkga`](@ref), and write them using
[`write_configuration`](@ref).

Fields
------
$(FIELDS)
"""
mutable struct BrkgaParams
    ########################################
    # BRKGA Hyper-parameters
    ########################################
    """
    Number of elements in the population [> 0].
    """
    population_size::Int64

    """
    Percentage of individuals to become the elite set (0, 1].
    """
    elite_percentage::Float64

    """
    Percentage of mutants to be inserted in the population
    """
    mutants_percentage::Float64

    """
    Number of elite parents for mating [> 0].
    """
    num_elite_parents::Int64

    """
    Number of total parents for mating [> 0].
    """
    total_parents::Int64

    """
    Type of bias that will be used.
    """
    bias_type::BiasFunction

    """
    Number of independent parallel populations.
    """
    num_independent_populations::Int64

    ########################################
    # Path Relinking parameters
    ########################################

    """
    Number of pairs of chromosomes to be tested to path relinking.
    """
    pr_number_pairs::Int64

    """
    Mininum distance between chromosomes selected to path-relinking.
    """
    pr_minimum_distance::Float64

    """
    Path relinking type.
    """
    pr_type::PathRelinkingType

    """
    Individual selection to path-relinking.
    """
    pr_selection::PathRelinkingSelection

    """
    Defines the block size based on the size of the population.
    """
    alpha_block_size::Float64

    """
    Percentage / path size to be computed. Value in (0, 1].
    """
    pr_percentage::Float64

    """
    Initialization constructor.
    """
    BrkgaParams() = new(0, 0, 0, 0, 0, CONSTANT, 0, 0, 0.0, DIRECT,
                        BESTSOLUTION, 0.0, 0.0)
end

################################################################################

"""
    mutable struct ExternalControlParams

Represents additional control parameters that can be used outside this
framework. You can load these parameters from a configuration file using
[`load_configuration`](@ref) and [`build_brkga`](@ref), and write them using
[`write_configuration`](@ref).

Fields
------
$(FIELDS)
"""
mutable struct ExternalControlParams
    """
    Interval at which elite chromosomes are exchanged
    (0 means no exchange) [> 0].
    """
    exchange_interval::Int64

    """
    Number of elite chromosomes exchanged from each population [> 0].
    """
    num_exchange_indivuduals::Int64

    """
    Interval at which the populations are reset (0 means no reset) [> 0].
    """
    reset_interval::Int64

    """
    Initialization constructor.
    """
    ExternalControlParams(
        exchange_interval::Int64 = 0,
        num_exchange_indivuduals::Int64 = 0,
        reset_interval::Int64 = 0) = new(exchange_interval,
                                         num_exchange_indivuduals,
                                         reset_interval)
end

################################################################################

"""
    mutable struct BrkgaData

Represents the internal state of the BRKGA-MP-IPR algorithm.

This structure has no direct constructor and must be built using
[`build_brkga`](@ref) functions. You can create multiple `BrkgaData`
representing different states of the algorithm, and use them independently.

!!! warning
    This structure is **NOT INTENDED** to be used outside BRKGA functions.
    Ad hoc changes may lead to inadvertent results.

Fields
------
$(FIELDS)
"""
mutable struct BrkgaData
    ########################################
    # Hyper-parameters
    ########################################

    """
    Optimization sense (minimization = 0, maximization = 1).
    """
    opt_sense::Sense

    """
    Number of genes in the chromosome [> 0].
    """
    chromosome_size::Int64

    """
    BRKGA parameters for evolution.
    """
    params::BrkgaParams

    """
    Number of elite items in the population [> 0].
    """
    elite_size::Int64

    """
    Number of mutants introduced at each generation into the population [> 0].
    """
    num_mutants::Int64

    """
    If false, no evolution is performed but only chromosome decoding.
    Very useful to emulate a multi-start algorithm.
    """
    evolutionary_mechanism_on::Bool

    ########################################
    # Internal data
    ########################################

    """
    *(Internal data)*
    The problem instance used by the `decode!` function to map a chromosome
    to a problem solution. Since `decode!` should not change this data,
    this attribute can be considered constant.
    """
    problem_instance::AbstractInstance

    """
    *(Internal data)*
    This is the **main decode function** called during the evolutionary
    process and in the path relink. It **must have** the following signature:

    ```julia
    decode!(chromosome::Array{Float64, 1},
            problem_instance::AbstractInstance,
            rewrite::Bool = true)::Float64
    ```

    Note that if `rewrite == false`, `decode!` must not change `chromosome`.
    IPR routines requires `decode!` to not change `chromosome`.
    """
    # TODO (ceandrade): the current Julia version (1.0) doesn't support
    # strong typing function signatures, as defined above. When such cabability
    # be available in Julia, the definition below must be changed for a
    # strong typed version.
    decode!::Function

    """
    *(Internal data)*
    The internal random generator number. DON'T USE FOR ANYTHING OUTSIDE.
    If you need a RNG, create a new generator.
    """
    rng::Random.MersenneTwister

    """
    *(Internal data)*
    Previous population.
    """
    previous::Array{Population, 1}

    """
    *(Internal data)*
    Current population.
    """
    current::Array{Population, 1}

    """
    *(Internal data)*
    A unary non-increasing function such that
    `bias_function(::Int64 > 0)::Float64`
    """
    bias_function::Function

    """
    *(Internal data)*
    Holds the sum of the results of each raking given a bias function.
    This value is needed to normalization.
    """
    total_bias_weight::Float64

    """
    *(Internal data)*
    Used to shuffled individual/chromosome indices during the mate.
    """
    shuffled_individuals::Array{Int64, 1}

    """
    *(Internal data)*
    Defines the order of parents during the mating.
    """
    parents_ordered::Array{Tuple{Float64, Int64}, 1}

    """
    *(Internal data)*
    Indicates if the algorithm was proper initialized.
    """
    initialized::Bool

    """
    *(Internal data)*
    Indicates if the algorithm have been reset.
    """
    reset_phase::Bool
end
