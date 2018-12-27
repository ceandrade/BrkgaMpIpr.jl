################################################################################
# types.jl: Definitions of internal data structures and external API.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Dec 26, 2018 by ceandrade
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

Specifies the objective function as minimization (0) or maximization (1).
"""
@enum Sense begin
   MINIMIZE = 0
   MAXIMIZE = 1
end

################################################################################

"""
    @enum BiasFunction

Specifies a bias function when choosing parents to mating. This function
substitutes the `rhoe` parameter from the original BRKGA. For a given
rank r, we have the following functions:

    - CONSTANT: 1 / number of parents for mating (all individuals have the
                same probability)

    - CUBIC: r^-3

    - EXPONENTIAL: ϵ^-r

    - LINEAR: 1/r

    - LOGINVERSE: 1 / log(r + 1)

    - QUADRATIC: r^-2
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

    - DIRECT: changes each key for the correspondent in the other chromosome.

    - PERMUTATION: switches the order of a key for that in the other chromosome.
"""
@enum PathRelinkingType begin
    DIRECT = 0
    PERMUTATION = 1
end

################################################################################

"""
    @enum PathRelinkingSelection

Specifies which individuals used to build the path:

    - BESTSOLUTION: selects, in the order, the best solution of each population.

    - RANDOMELITE: chooses uniformly random solutions from the elite sets.
"""
@enum PathRelinkingSelection begin
    BESTSOLUTION = 0
    RANDOMELITE = 1
end

################################################################################

"""
    @enum PathRelinkingResult

Specifies the result type/status of path relink procedure:

    - TOO_HOMOGENEOUS: the chromosomes among the populations are too homogeneous
                       and the path relink will not generate improveded
                       solutions.

    - NO_IMPROVEMENT: path relink was done but no improveded solution was found.

    - ELITE_IMPROVEMENT: an improved solution among the elite set was found,
                         but the best solution was not improved.

    - BEST_IMPROVEMENT: the best solution was improved.
"""
@enum PathRelinkingResult begin
    TOO_HOMOGENEOUS = 0
    NO_IMPROVEMENT = 1
    ELITE_IMPROVEMENT = 2
    BEST_IMPROVEMENT = 3
end

################################################################################
# Data structures
################################################################################

"""
    const empty_function() = nothing

Represent an empty function to be used as flag during data and bias function
setups.
"""
const empty_function() = nothing

################################################################################

"""
    abstract type AbstractInstance

The required interface for external data to be provided to the decoder. The
problem definitions and data **must be** a subtype of AbstractInstance. For
example,

```julia
mutable struct TSPInstance <: AbstractInstance
    num_cities::Int64
    distances::Array{Float64,2}
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
    mutable struct BrkgaData

Represents the internal state of the BRKGA-MP-IPR algorithm. This structure is
**NOT intended** to be used outside BRKGA functions. Ad hoc changes may lead to
inadvertent results. However, some fields of interest are documented using
docstring. Internal only fields have regular comments as documentation.
"""
mutable struct BrkgaData
    ########################################
    # BRKGA Hyper-parameter
    ########################################

    """
    (BRKGA Hyper-parameter)
    Optimization sense (minimization = 0, maximization = 1).
    """
    opt_sense::Sense

    """
    (BRKGA Hyper-parameter)
    Number of genes in the chromosome [> 0].
    """
    chromosome_size::Int64

    """
    (BRKGA Hyper-parameter)
    Number of elements in the population [> 0].
    """
    population_size::Int64

    """
    (BRKGA Hyper-parameter)
    Number of elite items in the population [> 0].
    """
    elite_size::Int64

    """
    (BRKGA Hyper-parameter)
    Number of mutants introduced at each generation into the population [> 0].
    """
    num_mutants::Int64

    """
    (BRKGA Hyper-parameter)
    Number of elite parents for mating [> 0].
    """
    num_elite_parents::Int64

    """
    (BRKGA Hyper-parameter)
    Number of total parents for mating [> 0].
    """
    total_parents::Int64

    """
    (BRKGA Hyper-parameter)
    Type of bias that will be used.
    """
    bias::BiasFunction

    """
    (BRKGA Hyper-parameter)
    Number of independent parallel populations.
    """
    num_independent_populations::Int64

    """
    (BRKGA Hyper-parameter)
    If false, no evolution is performed but only chromosome decoding.
    Very useful to emulate a multi-start algorithm.
    """
    evolutionary_mechanism_on::Bool

    ########################################
    # Path Relinking parameters
    ########################################

    # TODO (ceandrade): list the path relink parameters here.

    ########################################
    # Internal data
    ########################################

    """
    (Internal data)
    The problem instance used by the `decode!` function to map a chromosome
    to a problem solution. Since `decode!` should not change this data,
    this attribute can be considered constant.
    """
    problem_instance::AbstractInstance

    """
    (Internal data)
    This is the **main decode function** called during the evolutionary
    process and in the path relink. It **must have** the following signature:

        decode!(chromosome::Array{Float64, 1},
                problem_instance::AbstractInstance,
                rewrite::Bool = true)::Float64

    Note that if `rewrite == true`, `decode!` can change `chromosome`.
    """
    # TODO (ceandrade): the current Julia version (0.6) doesn't support
    # strong typing function signatures, as defined above. When such cabability
    # be available in Julia, the definition below must be changed for a
    # strong typed version.
    decode!::Function

    """
    (Internal data)
    The internal random generator number. DON'T USE FOR ANYTHING OUTSIDE.
    If you need a RNG, create a new generator.
    """
    rng::Random.MersenneTwister

    """
    (Internal data)
    Previous population.
    """
    previous::Array{Population, 1}

    """
    (Internal data)
    Current population.
    """
    current::Array{Population, 1}

    """
    (Internal data)
    A unary non-increasing function such that
    `bias_function(::Int64 > 0)::Float64`
    """
    bias_function::Function

    """
    (Internal data)
    Holds the sum of the results of each raking given a bias function.
    This value is needed to normalization.
    """
    total_bias_weight::Float64

    """
    (Internal data)
    Used to shuffled individual/chromosome indices during the mate.
    """
    shuffled_individuals::Array{Int64, 1}

    """
    (Internal data)
    Defines the order of parents during the mating.
    """
    parents_ordered::Array{Tuple{Float64, Int64}, 1}

    """
    (Internal data)
    Indicates if the algorithm was proper initialized.
    """
    initialized::Bool

    """
    (Internal data)
    Indicates if the algorithm have been reset.
    """
    reset_phase::Bool
end

################################################################################

"""
    struct ControlParams

Represents additional control parameters that can be used outside this
framework. They can be loaded from a configuration file by `init()`.
"""
struct ExternalControlParams
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
