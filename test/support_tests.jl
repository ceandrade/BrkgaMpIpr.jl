################################################################################
# support_tests.jl: unit tests for support routines of BrkgaMpIpr.
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

@testset "get_best_fitness()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["seed"]] = 12323
    param_values[param_index["chr_size"]] = CHROMOSOME_SIZE
    param_values[param_index["brkga_params"]].population_size = 5000
    param_values[param_index["brkga_params"]].num_independent_populations = 10

    brkga_data = build_brkga(param_values...)
    params = brkga_data.params

    # Not initialized
    @test_throws ErrorException get_best_fitness(brkga_data)

    initialize!(brkga_data)

    ########################
    # Test for maximization
    ########################

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    num_individuals = params.population_size *
                      params.num_independent_populations

    best_value = -Inf
    for i in 1:num_individuals
        local_chr = rand(local_rng, param_values[param_index["chr_size"]])
        best_value = max(best_value, sum_decode!(local_chr, instance))
    end

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    @test get_best_fitness(brkga_data) ≈ best_value

    ########################
    # Test for minimization
    ########################

    param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params
    initialize!(brkga_data)

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    num_individuals = params.population_size *
                      params.num_independent_populations

    best_value = Inf
    for i in 1:num_individuals
        local_chr = rand(local_rng, param_values[param_index["chr_size"]])
        best_value = min(best_value, sum_decode!(local_chr, instance))
    end

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    @test get_best_fitness(brkga_data) ≈ best_value
end

################################################################################

@testset "get_best_chromosome()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["seed"]] = 12323
    param_values[param_index["chr_size"]] = CHROMOSOME_SIZE
    param_values[param_index["brkga_params"]].population_size = 5000
    param_values[param_index["brkga_params"]].num_independent_populations = 10

    brkga_data = build_brkga(param_values...)
    params = brkga_data.params

    # Not initialized
    @test_throws ErrorException get_best_chromosome(brkga_data)

    initialize!(brkga_data)

    ########################
    # Test for maximization
    ########################

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    num_individuals = params.population_size *
                      params.num_independent_populations

    best_value = -Inf
    best_chr = Array{Float64, 1}()
    for i in 1:num_individuals
        local_chr = rand(local_rng, param_values[param_index["chr_size"]])
        value = sum_decode!(local_chr, instance)

        if value > best_value
            best_value = value
            best_chr = local_chr
        end
    end

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    @test get_best_chromosome(brkga_data) ≈ best_chr

    ########################
    # Test for minimization
    ########################

    param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)
    params = brkga_data.params
    initialize!(brkga_data)

    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)

    best_value = Inf
    best_chr = Array{Float64, 1}()
    for i in 1:num_individuals
        local_chr = rand(local_rng, param_values[param_index["chr_size"]])
        value = sum_decode!(local_chr, instance)
        if value < best_value
            best_value = value
            best_chr = local_chr
       end
    end

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    @test get_best_chromosome(brkga_data) ≈ best_chr
end

################################################################################

@testset "get_chromosome()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_independent_populations = 3
    brkga_data = build_brkga(param_values...)

    # Not initialized
    @test_throws ErrorException get_chromosome(brkga_data, 1, 1)

    initialize!(brkga_data)

    # Test invalid population indices.
    @test_throws ArgumentError get_chromosome(brkga_data, 0, 1)
    @test_throws ArgumentError get_chromosome(brkga_data,
        brkga_data.params.num_independent_populations + 1, 1)

    # Test invalid chrmosome indices.
    @test_throws ArgumentError get_chromosome(brkga_data, 1, 0)
    @test_throws ArgumentError get_chromosome(brkga_data, 1,
        brkga_data.params.population_size + 1)

    # TODO (ceandrade): this test is not correct. Please, replicate the Python
    # tests here.

    idx, pop = 1, 1
    actual_chr = brkga_data.current[pop]
    copy_chr = get_chromosome(brkga_data, pop, idx)
    @test copy_chr !== actual_chr.chromosomes[actual_chr.fitness[idx][2]]
    @test copy_chr ≈ actual_chr.chromosomes[actual_chr.fitness[idx][2]]

    idx, pop = 2, 2
    actual_chr = brkga_data.current[pop]
    copy_chr = get_chromosome(brkga_data, pop, idx)
    @test copy_chr !== actual_chr.chromosomes[actual_chr.fitness[idx][2]]
    @test copy_chr ≈ actual_chr.chromosomes[actual_chr.fitness[idx][2]]

    idx = brkga_data.params.population_size
    pop = brkga_data.params.num_independent_populations
    actual_chr = brkga_data.current[pop]
    copy_chr = get_chromosome(brkga_data, pop, idx)
    @test copy_chr !== actual_chr.chromosomes[actual_chr.fitness[idx][2]]
    @test copy_chr ≈ actual_chr.chromosomes[actual_chr.fitness[idx][2]]
end

################################################################################

@testset "get_current_population()" begin
    param_values = deepcopy(default_param_values)
    param_values[param_index["brkga_params"]].num_independent_populations = 3
    brkga_data = build_brkga(param_values...)

    # Not initialized
    @test_throws ErrorException get_current_population(brkga_data, 1)

    initialize!(brkga_data)

    @test_throws ArgumentError get_current_population(brkga_data, 0)
    @test_throws ArgumentError get_current_population(brkga_data,
        brkga_data.params.num_independent_populations + 1)

    for i in 1:brkga_data.params.num_independent_populations
        @test get_current_population(brkga_data, i) === brkga_data.current[i]
    end
end
