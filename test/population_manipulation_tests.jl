################################################################################
# population_manipulation_tests.jl: unit tests for support routines of
#                                   BrkgaMpIpr.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Feb 28, 2019 by ceandrade
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

@testset "reset!()" begin
    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)

    @test_throws ErrorException reset!(brkga_data)

    initialize!(brkga_data)

    # Create a local RNG and advance it until the same state as the internal
    # BrkgaData RNG after initialization.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    skip = brkga_data.params.num_independent_populations *
           brkga_data.params.population_size * brkga_data.chromosome_size

    # Assert the both generators are in the same state.
    rand(local_rng, 1000 + skip)
    @assert rand(brkga_data.rng) == rand(local_rng)

    # Create a local chromosome and applied the decoder on it.
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    sum_decode!(local_chr, instance)

    # Reset and test the first individual.
    reset!(brkga_data)
    @test brkga_data.current[1].chromosomes[1] == local_chr

    # After reset, the reset phase flag should be deactivated.
    @test brkga_data.reset_phase == false;
end

################################################################################

@testset "exchange_elite!()" begin
    ########################
    # Exceptions
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].population_size = 10
    param_values[param_index["brkga_params"]].num_independent_populations = 2
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params

    # Not initialized
    @test_throws ErrorException exchange_elite!(brkga_data, 1)

    initialize!(brkga_data)

    # Wrong number of individuals to exchange.
    @test_throws ArgumentError exchange_elite!(brkga_data, 0)
    @test_throws ArgumentError exchange_elite!(brkga_data, -10)
    @test_throws ArgumentError exchange_elite!(brkga_data, params.population_size)
    @test_throws ArgumentError exchange_elite!(brkga_data, params.population_size + 10)

    # More exchanges than number of chromosomes.
    num_immigrants = cld(params.population_size,
                         params.num_independent_populations - 1)

    @test_throws ArgumentError exchange_elite!(brkga_data, num_immigrants)
    @test_throws ArgumentError exchange_elite!(brkga_data, num_immigrants + 1)

    ########################
    # Single population, no exchange
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_independent_populations = 1
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    local_chromosomes = deepcopy(brkga_data.current[1].chromosomes)
    local_fitness = deepcopy(brkga_data.current[1].fitness)

    exchange_elite!(brkga_data, 1)

    @test local_chromosomes == brkga_data.current[1].chromosomes
    @test local_fitness == brkga_data.current[1].fitness

    ########################
    # Two populations, maximization
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["brkga_params"]].population_size = 10
    param_values[param_index["brkga_params"]].num_independent_populations = 2
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    if brkga_data.current[1].fitness[1][1] > brkga_data.current[2].fitness[1][1]
        pop = 1
    else
        pop = 2
    end

    value_idx = brkga_data.current[pop].fitness[1]
    best_value = value_idx[1]
    best_chr = copy(brkga_data.current[pop].chromosomes[value_idx[2]])

    exchange_elite!(brkga_data, 1)

    value_idx = brkga_data.current[1].fitness[1]
    @test value_idx[1] == best_value
    @test brkga_data.current[1].chromosomes[value_idx[2]] == best_chr

    value_idx = brkga_data.current[2].fitness[1]
    @test value_idx[1] == best_value
    @test brkga_data.current[2].chromosomes[value_idx[2]] == best_chr

    ########################
    # Five populations, minimization
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["instance"]] = Instance(1000)
    param_values[param_index["chr_size"]] = 1000
    param_values[param_index["opt_sense"]] = MINIMIZE
    param_values[param_index["brkga_params"]].population_size = 1000
    param_values[param_index["brkga_params"]].num_independent_populations = 5
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    best_value = brkga_data.current[1].fitness[end][1]
    pop = 1
    for i in 1:brkga_data.params.num_independent_populations
        if best_value > brkga_data.current[i].fitness[1][1]
            best_value = brkga_data.current[i].fitness[1][1]
            pop = i
        end
    end

    value_idx = brkga_data.current[pop].fitness[1]
    best_value = value_idx[1]
    best_chr = copy(brkga_data.current[pop].chromosomes[value_idx[2]])

    exchange_elite!(brkga_data, 1)

    for i in 1:brkga_data.params.num_independent_populations
        value_idx = brkga_data.current[i].fitness[1]
        @test value_idx[1] == best_value
        @test brkga_data.current[i].chromosomes[value_idx[2]] == best_chr
    end
end

################################################################################

@testset "inject_chromosome!()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["brkga_params"]].num_independent_populations = 3
    brkga_data = build_brkga(param_values...)

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    local_chr = rand(local_rng, brkga_data.chromosome_size)

    # Not initialized
    @test_throws ErrorException inject_chromosome!(brkga_data, local_chr, 1,
                                                   1, 0.0)

    initialize!(brkga_data)

    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr, 0, 1, 0.0)
    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr,
        brkga_data.params.num_independent_populations + 1, 1, 0.0)

    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr, 1, 0, 0.0)
    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr, 1,
        brkga_data.params.population_size + 1, 0.0)

    local_chr = rand(local_rng, brkga_data.chromosome_size - 1)
    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr, 1, 1, 0.0)

    local_chr = rand(local_rng, brkga_data.chromosome_size + 1)
    @test_throws ArgumentError inject_chromosome!(brkga_data, local_chr, 1, 1, 0.0)

    # This should create the best solution after decoding with value 100.0
    max = get_best_fitness(brkga_data)
    local_chr = fill(max, brkga_data.chromosome_size)
    local_chr -= brkga_data.problem_instance.data

    # Insert a chromosome in the last position.
    inject_chromosome!(brkga_data, local_chr, 1,
                       brkga_data.params.population_size)
    @test get_best_fitness(brkga_data) ≈ 100.0

    # Insert one a little better in the middle.
    inject_chromosome!(brkga_data, local_chr, 1,
                       brkga_data.params.population_size ÷ 2, 101.0)
    @test get_best_fitness(brkga_data) ≈ 101.0

    # Insert a bad one in the middle.
    inject_chromosome!(brkga_data, local_chr, 1,
                       brkga_data.params.population_size ÷ 2, -10.0)
    @test get_best_fitness(brkga_data) ≈ 101.0

    # Since is maximization, the bad one must be in the last position.
    @test get_best_fitness(brkga_data) != -10.0
    @test brkga_data.current[1].fitness[end][1] ≈ -10.0

    # Population 2 must be intact.
    @test brkga_data.current[2].fitness[1][1] != 101.0
    @test brkga_data.current[2].fitness[end][1] != -10.0
end

################################################################################

@testset "shake!()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["brkga_params"]].population_size = 5
    param_values[param_index["brkga_params"]].elite_percentage = 0.4
    param_values[param_index["brkga_params"]].num_independent_populations = 1
    param_values[param_index["instance"]] = Instance(10)
    default_param_values[param_index["decode!"]] = sum_decode!

    brkga_data = build_brkga(param_values...)

    @test_throws ErrorException shake!(brkga_data, 1, CHANGE, 1)

    initialize!(brkga_data)

    @assert brkga_data.elite_size == 2

    @test_throws ArgumentError shake!(brkga_data, -1, CHANGE, 1)
    @test_throws ArgumentError shake!(brkga_data, 1, CHANGE, -1)

    # TODO (ceandrade): implement the full tests.
end
