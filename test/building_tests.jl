################################################################################
# building_tests.jl: unit tests for building routines of BrkgaMpIpr.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Jan 04, 2019 by ceandrade
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

    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)

    @test brkga_data.elite_size == 3
    @test brkga_data.num_mutants == 1

    brkga_params = param_values[param_index["brkga_params"]]
    @test length(brkga_data.shuffled_individuals) == brkga_params.population_size
    @test length(brkga_data.parents_ordered) == brkga_params.total_parents

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    # Same warm up that in build_brkga().
    rand(local_rng, 1000)
    @test rand(brkga_data.rng) == rand(local_rng)

    ########################
    # Test multi-start building.
    ########################

    param_values[param_index["evolutionary_mechanism_on"]] = false
    param_values[param_index["brkga_params"]].population_size = 10
    brkga_data = build_brkga(param_values...)
    @test brkga_data.elite_size == 1
    @test brkga_data.num_mutants == 9

    ########################
    # Test bias functions.
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].bias_type = LOGINVERSE
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.4426950408889634
    @test brkga_data.bias_function(2) ≈ 0.9102392266268375
    @test brkga_data.bias_function(3) ≈ 0.7213475204444817

    param_values[param_index["brkga_params"]].bias_type = LINEAR
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.5
    @test brkga_data.bias_function(3) ≈ 0.3333333333333333

    param_values[param_index["brkga_params"]].bias_type = QUADRATIC
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.25
    @test brkga_data.bias_function(3) ≈ 0.1111111111111111

    param_values[param_index["brkga_params"]].bias_type = CUBIC
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 1.0
    @test brkga_data.bias_function(2) ≈ 0.125
    @test brkga_data.bias_function(3) ≈ 0.037037037037037035

    param_values[param_index["brkga_params"]].bias_type = EXPONENTIAL
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 0.36787944117144233
    @test brkga_data.bias_function(2) ≈ 0.1353352832366127
    @test brkga_data.bias_function(3) ≈ 0.049787068367863944

    param_values[param_index["brkga_params"]].bias_type = CONSTANT
    brkga_data = build_brkga(param_values...)
    @test brkga_data.bias_function(1) ≈ 0.5
    @test brkga_data.bias_function(2) ≈ 0.5
    @test brkga_data.bias_function(3) ≈ 0.5

    ########################
    # Test exceptions.
    ########################

    # Chromosome size
    param_values = deepcopy(default_param_values)
    param_values[param_index["chr_size"]] = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["chr_size"]] = -10
    @test_throws ArgumentError build_brkga(param_values...)

    # Population size
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].population_size = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].population_size = -10
    @test_throws ArgumentError build_brkga(param_values...)

    # Elite size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].elite_percentage = 0.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].elite_percentage = -1.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].elite_percentage = 0.3
    param_values[param_index["brkga_params"]].population_size = 2
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].elite_percentage = 1.1
    param_values[param_index["brkga_params"]].population_size = 10
    @test_throws ArgumentError build_brkga(param_values...)

    # Mutant size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].mutants_percentage = -1.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].mutants_percentage = 1.0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].mutants_percentage = 1.1
    @test_throws ArgumentError build_brkga(param_values...)

    # Elite + Mutant size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].elite_percentage = 0.6
    param_values[param_index["brkga_params"]].mutants_percentage = 0.6
    @test_throws ArgumentError build_brkga(param_values...)

    # Elite parents for mating.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_elite_parents = 0
    @test_throws ArgumentError build_brkga(param_values...)

    param_values[param_index["brkga_params"]].num_elite_parents = 2
    param_values[param_index["brkga_params"]].total_parents = 2
    @test_throws ArgumentError build_brkga(param_values...)

    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_elite_parents = 1 +
            ceil(Int64, brkga_params.population_size *
                        brkga_params.elite_percentage)
    brkga_params.total_parents = 1 + brkga_params.num_elite_parents
    @test_throws ArgumentError build_brkga(param_values...)

    # Number of independent populations.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_independent_populations = 0
    @test_throws ArgumentError build_brkga(param_values...)

    # alpha_block_size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].alpha_block_size = 0.0
    @test_throws ArgumentError build_brkga(param_values...)

    # percentage / path size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].pr_percentage = 0.0
    @test_throws ArgumentError build_brkga(param_values...)

    # percentage / path size.
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].pr_percentage = 1.001
    @test_throws ArgumentError build_brkga(param_values...)
end

################################################################################

@testset "Config. file build_brkga()" begin
    config_path = joinpath(@__DIR__, "configuration_files")

    local_param_values = [
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

    @test brkga_data.elite_size == 150
    @test brkga_data.num_mutants == 75
    @test brkga_data.params.population_size == 500
    @test brkga_data.params.num_elite_parents == 2
    @test brkga_data.params.total_parents == 3
    @test brkga_data.params.bias_type == LOGINVERSE
    @test brkga_data.params.num_independent_populations == 3
    @test brkga_data.params.pr_number_pairs == 0
    @test brkga_data.params.pr_minimum_distance == 0.15
    @test brkga_data.params.pr_type == PERMUTATION
    @test brkga_data.params.pr_selection == RANDOMELITE
    @test brkga_data.params.alpha_block_size == 1.0
    @test brkga_data.params.pr_percentage == 1.0
    @test external_params.exchange_interval == 200
    @test external_params.num_exchange_indivuduals == 2
    @test external_params.reset_interval == 600
end

################################################################################

@testset "initialize!()" begin
    ########################
    # Test with custom function
    # loaded from configuration file
    ########################

    config_path = joinpath(@__DIR__, "configuration_files")
    local_param_values = [
        default_param_values[param_index["instance"]],
        default_param_values[param_index["decode!"]],
        default_param_values[param_index["opt_sense"]],
        default_param_values[param_index["seed"]],
        default_param_values[param_index["chr_size"]]
    ]

    brkga_data, external_params =
        build_brkga(local_param_values...,
                    joinpath(config_path, "custom_bias_function.conf"))

    # Custom function is not defined.
    @test_throws ErrorException initialize!(brkga_data)

    ########################
    # Test without warmstart
    ########################
    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params

    initialize!(brkga_data)

    for i in 1:params.num_independent_populations
        @test length(brkga_data.current) == params.num_independent_populations
        @test length(brkga_data.current[i].chromosomes) == params.population_size
        @test length(brkga_data.current[i].fitness) == params.population_size

        @test length(brkga_data.previous) == params.num_independent_populations
        @test length(brkga_data.previous[i].chromosomes) == params.population_size
        @test length(brkga_data.previous[i].fitness) == params.population_size

        @test brkga_data.current[i].chromosomes == brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].chromosomes !== brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].fitness == brkga_data.previous[i].fitness
        @test brkga_data.current[i].fitness !== brkga_data.previous[i].fitness

        correct_order = true
        for j in 2:length(brkga_data.current[i].fitness)
            correct_order &= brkga_data.current[i].fitness[j-1] >=
                             brkga_data.current[i].fitness[j]
        end
        @test correct_order
    end

    @test brkga_data.initialized == true
    @test brkga_data.reset_phase == false

    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params
    initialize!(brkga_data)

    for i in 1:params.num_independent_populations
        correct_order = true
        for j in 2:length(brkga_data.current[i].fitness)
            correct_order &= brkga_data.current[i].fitness[j-1] <=
                             brkga_data.current[i].fitness[j]
        end
        @test correct_order
    end

    ########################
    # Test with warmstart
    ########################

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    chromosomes = [
        rand(local_rng, param_values[param_index["chr_size"]]),
        rand(local_rng, param_values[param_index["chr_size"]]),
        rand(local_rng, param_values[param_index["chr_size"]])
    ]
    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params
    set_initial_population!(brkga_data, chromosomes)

    initialize!(brkga_data)

    for i in 1:params.num_independent_populations
        @test length(brkga_data.current) == params.num_independent_populations
        @test length(brkga_data.current[i].chromosomes) == params.population_size
        @test length(brkga_data.current[i].fitness) == params.population_size

        @test length(brkga_data.previous) == params.num_independent_populations
        @test length(brkga_data.previous[i].chromosomes) == params.population_size
        @test length(brkga_data.previous[i].fitness) == params.population_size

        @test brkga_data.current[i].chromosomes == brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].chromosomes !== brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].fitness == brkga_data.previous[i].fitness
        @test brkga_data.current[i].fitness !== brkga_data.previous[i].fitness
    end

    sum_decode!(chromosomes[1], instance)
    @test brkga_data.current[1].chromosomes[1] == chromosomes[1]

    # Create a local chromosome and applied the decoder on it.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    sum_decode!(local_chr, instance)

    # 4th chromosome must be the 1st generated due to the warmstart.
    @test brkga_data.current[1].chromosomes[4] == local_chr

    ########################
    # Test reset phase
    ########################

    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params
    initialize!(brkga_data)

    # Create a local RNG and advance it until the same state as the internal
    # BrkgaData RNG after initialization.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    skip = params.num_independent_populations *
           params.population_size * brkga_data.chromosome_size
    rand(local_rng, 1000 + skip)

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    brkga_data.reset_phase = true
    initialize!(brkga_data)

    # Create a local chromosome and applied the decoder on it.
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    sum_decode!(local_chr, instance)

    @test brkga_data.current[1].chromosomes[1] == local_chr
end

###############################################################################

@testset "set_bias_custom_function!()" begin
    param_values = deepcopy(default_param_values)

    param_values[param_index["brkga_params"]].population_size = 100
    param_values[param_index["brkga_params"]].total_parents = 10

    brkga_data = build_brkga(param_values...)

    # After build, brkga_params function is never CUSTOM
    @test brkga_data.params.bias_type != CUSTOM

    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> x)
    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> x + 1)
    @test_throws ArgumentError set_bias_custom_function!(brkga_data, x -> log1p(x))

    set_bias_custom_function!(brkga_data, x -> 1.0 / log1p(x))
    @test brkga_data.total_bias_weight ≈ 6.554970525044798

    # After 2nd call to set_bias_custom_function, brkga_params function is always CUSTOM
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 1.0 / x)
    @test brkga_data.total_bias_weight ≈ 2.9289682539682538
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> x ^ -2.0)
    @test brkga_data.total_bias_weight ≈ 1.5497677311665408
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> x ^ -3.0)
    @test brkga_data.total_bias_weight ≈ 1.197531985674193
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> exp(-x))
    @test brkga_data.total_bias_weight ≈ 0.5819502851677112
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 1.0 / brkga_data.params.total_parents)
    @test brkga_data.total_bias_weight ≈ 0.9999999999999999
    @test brkga_data.params.bias_type == CUSTOM

    set_bias_custom_function!(brkga_data, x -> 0.6325 / sqrt(x))
    @test brkga_data.total_bias_weight ≈ 3.175781171302612
    @test brkga_data.params.bias_type == CUSTOM

    #############################################
    # Constant functions test for standard BRKGA
    #############################################

    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_elite_parents = 1
    param_values[param_index["brkga_params"]].total_parents = 2
    brkga_data = build_brkga(param_values...)

    rho = 0.5
    set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)
    @test brkga_data.total_bias_weight ≈ 1.0

    rho = 0.75
    set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)
    @test brkga_data.total_bias_weight ≈ 1.0

    rho = 0.90
    set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)
    @test brkga_data.total_bias_weight ≈ 1.0
end

################################################################################

@testset "set_initial_population!()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["chr_size"]] = 3
    param_values[param_index["brkga_params"]].num_independent_populations = 2
    brkga_data = build_brkga(param_values...)
    local_rng = MersenneTwister(param_values[param_index["seed"]])

    chromosomes = Array{Array{Float64, 1}, 1}(undef,
        param_values[param_index["brkga_params"]].population_size + 1
    )
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 1)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] - 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]] + 1)
    @test_throws ArgumentError set_initial_population!(brkga_data, chromosomes)

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 1)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    set_initial_population!(brkga_data, chromosomes)

    @test length(brkga_data.current[1].chromosomes) == length(chromosomes)
    @test brkga_data.current[1].chromosomes == chromosomes
    @test brkga_data.current[1].chromosomes !== chromosomes

    chromosomes[1] = [0.1111, 0.2222, 0.3333]
    @test brkga_data.current[1].chromosomes != chromosomes

    chromosomes = Array{Array{Float64, 1}, 1}(undef, 3)
    chromosomes[1] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[2] = rand(local_rng, param_values[param_index["chr_size"]])
    chromosomes[3] = rand(local_rng, param_values[param_index["chr_size"]])
    set_initial_population!(brkga_data, chromosomes)
    @test length(brkga_data.current[1].chromosomes) == length(chromosomes)
    @test brkga_data.current[1].chromosomes == chromosomes
    @test brkga_data.current[1].chromosomes !== chromosomes
end
