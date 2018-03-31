################################################################################
# building_tests.jl: unit tests for building routines of BrkgaMpIpr.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Mar 31, 2018 by ceandrade
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

@testset "Detailed build_brkga()" begin
    ########################
    # Test regular/correct building.
    ########################

    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)
    @test brkga_data.elite_size == 3
    @test brkga_data.num_mutants == 1
    # @test length(brkga_data.previous) == param_values[param_index["num_independent_populations"]]
    # @test length(brkga_data.current) == param_values[param_index["num_independent_populations"]]
    @test length(brkga_data.shuffled_individuals) == param_values[param_index["pop_size"]]
    @test length(brkga_data.parents_ordered) == param_values[param_index["pop_size"]]

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    # Same warm up that in build_brkga().
    rand(local_rng, 1000)
    @test rand(brkga_data.rng) == rand(local_rng)

    ########################
    # Test multi-start building.
    ########################

    param_values[param_index["evolutionary_mechanism_on"]] = false
    param_values[param_index["pop_size"]] = 10
    brkga_data = build_brkga(param_values...)
    @test brkga_data.elite_size == 1
    @test brkga_data.num_mutants == 9

    param_values = copy(default_param_values)

    ########################
    # Test bias functions.
    ########################

    param_values[param_index["bias"]] = LOGINVERSE
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.4426950408889634
    @test brkga_data.bias_function(2) ≈ 0.9102392266268375
    @test brkga_data.bias_function(3) ≈ 0.7213475204444817

    param_values[param_index["bias"]] = LINEAR
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.5
    @test brkga_data.bias_function(3) ≈ 0.3333333333333333

    param_values[param_index["bias"]] = QUADRATIC
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.25
    @test brkga_data.bias_function(3) ≈ 0.1111111111111111

    param_values[param_index["bias"]] = CUBIC
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.125
    @test brkga_data.bias_function(3) ≈ 0.037037037037037035

    param_values[param_index["bias"]] = EXPONENTIAL
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 0.36787944117144233
    @test brkga_data.bias_function(2) ≈ 0.1353352832366127
    @test brkga_data.bias_function(3) ≈ 0.049787068367863944

    param_values[param_index["bias"]] = CONSTANT
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 0.5
    @test brkga_data.bias_function(2) ≈ 0.5
    @test brkga_data.bias_function(3) ≈ 0.5

    param_values = copy(default_param_values)

    ########################
    # Test exceptions.
    ########################

    # Chromosome size
    param_values[param_index["chr_size"]] = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["chr_size"]] = -10
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Population size
    param_values[param_index["pop_size"]] = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["pop_size"]] = -10
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Elite size.
    param_values[param_index["elite_percentage"]] = 0.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["elite_percentage"]] = -10.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["elite_percentage"]] = 0.3
    param_values[param_index["pop_size"]] = 2
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["elite_percentage"]] = 1.1
    param_values[param_index["pop_size"]] = 10
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Mutant size.
    param_values[param_index["mutants_percentage"]] = -10.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["mutants_percentage"]] = 1.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["mutants_percentage"]] = 1.1
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Elite + Mutant size.
    param_values[param_index["elite_percentage"]] = 0.6
    param_values[param_index["mutants_percentage"]] = 0.6
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Elite parents for mating.
    param_values[param_index["num_elite_parents"]] = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["num_elite_parents"]] = 2
    param_values[param_index["total_parents"]] = 2
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["num_elite_parents"]] = 1 +
               ceil(Int64, param_values[param_index["pop_size"]] *
                           param_values[param_index["elite_percentage"]])
    param_values[param_index["total_parents"]] = 1 +
               param_values[param_index["num_elite_parents"]]
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # Number of independent populations.
    param_values[param_index["num_independent_populations"]] = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values = copy(default_param_values)

    # TODO (ceandrade): check path relink params here
end

################################################################################

@testset "Config. file build_brkga()" begin
    config_path = joinpath(@__DIR__, "configuration_files")

    const local_param_values = [
        default_param_values[param_index["instance"]],
        default_param_values[param_index["decode!"]],
        default_param_values[param_index["opt_sense"]],
        default_param_values[param_index["seed"]],
        default_param_values[param_index["chr_size"]]
    ]

    ########################
    # Test regular/correct building.
    ########################

    brkga_data, external_params =
        build_brkga(local_param_values...,
                    joinpath(config_path, "regular.conf"))

    @test brkga_data.population_size == 500
    @test brkga_data.elite_size == 150
    @test brkga_data.num_mutants == 75
    @test brkga_data.num_elite_parents == 2
    @test brkga_data.total_parents == 3
    @test brkga_data.bias == LOGINVERSE
    @test brkga_data.num_independent_populations == 3
    # TODO (ceandrade): test path relink parameters

    @test external_params.exchange_interval == 200
    @test external_params.num_exchange_indivuduals == 2
    @test external_params.reset_interval == 600

    ########################
    # Test exceptions.
    ########################

    @test_throws LoadError build_brkga(local_param_values..., ".")
    @test_throws SystemError build_brkga(local_param_values..., "")
    @test_throws LoadError build_brkga(local_param_values...,
            joinpath(config_path, "missing_value.conf"))
    @test_throws LoadError build_brkga(local_param_values...,
            joinpath(config_path, "unknown_param.conf"))
    @test_throws LoadError build_brkga(local_param_values...,
                    joinpath(config_path, "wrong_type.conf"))
    @test_throws LoadError build_brkga(local_param_values...,
                    joinpath(config_path, "missing_param.conf"))
    @test_throws LoadError build_brkga(local_param_values...,
                    joinpath(config_path, "wrong_bias_function.conf"))
end

################################################################################

@testset "set_bias_custom_function!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["total_parents"]] = 10
    brkga_data = build_brkga(param_values...)

    # After build, bias function is never CUSTOM
    @test brkga_data.bias != CUSTOM

    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> x)
    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> x + 1)
    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> log1p(x))

    set_bias_custom_function!(brkga_data, x -> 1.0 / log1p(x))
    @test brkga_data.total_bias_weight ≈ 6.554970525044798

    # After 2nd call to set_bias_custom_function, bias function is always CUSTOM
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 1.0 / x)
    @test brkga_data.total_bias_weight ≈ 2.9289682539682538
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> x ^ -2.0)
    @test brkga_data.total_bias_weight ≈ 1.5497677311665408
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> x ^ -3.0)
    @test brkga_data.total_bias_weight ≈ 1.197531985674193
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> exp(-x))
    @test brkga_data.total_bias_weight ≈ 0.5819502851677112
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 1.0 / brkga_data.total_parents)
    @test brkga_data.total_bias_weight ≈ 0.9999999999999999
    @test brkga_data.bias == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 0.6325 / sqrt(x))
    @test brkga_data.total_bias_weight ≈ 3.175781171302612
    @test brkga_data.bias == CUSTOM
end

################################################################################

@testset "set_initial_population!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["chr_size"]] = 3
    param_values[param_index["num_independent_populations"]] = 2
    brkga_data = build_brkga(param_values...)
    local_rng = MersenneTwister(param_values[param_index["seed"]])

    chromosomes = Array{Array{Float64, 1}, 1}(
        param_values[param_index["pop_size"]] + 1
    )
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(1)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] - 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(1)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    set_initial_population!(brkga_data, chromosomes)

    @test length(brkga_data.current[1].chromosomes) == length(chromosomes)
    @test brkga_data.current[1].chromosomes == chromosomes
    @test brkga_data.current[1].chromosomes !== chromosomes

    chromosomes[1] = [0.1111, 0.2222, 0.3333]
    @test brkga_data.current[1].chromosomes != chromosomes

    chromosomes = Array{Array{Float64, 1}, 1}(3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    set_initial_population!(brkga_data, chromosomes)
    @test length(brkga_data.current[1].chromosomes) == length(chromosomes)
    @test brkga_data.current[1].chromosomes == chromosomes
    @test brkga_data.current[1].chromosomes !== chromosomes
end
