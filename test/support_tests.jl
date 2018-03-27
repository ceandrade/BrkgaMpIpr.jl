################################################################################
# support_tests.jl: unit tests for support routines of BrkgaMpIpr.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
# Last update: Mar 27, 2018 by ceandrade
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

@testset "initialize!()" begin
    ########################
    # Test without warmstart
    ########################
    param_values = copy(default_param_values)
    param_values[param_index["opt_sense"]] = MAXIMIZE
    brkga_data = build_brkga(param_values...)

    initialize!(brkga_data)

    for i = 1:brkga_data.num_independent_populations
        @test length(brkga_data.current) == param_values[param_index["num_independent_populations"]]
        @test length(brkga_data.current[i].chromosomes) == param_values[param_index["pop_size"]]
        @test length(brkga_data.current[i].fitness) == param_values[param_index["pop_size"]]

        @test length(brkga_data.previous) == param_values[param_index["num_independent_populations"]]
        @test length(brkga_data.previous[i].chromosomes) == param_values[param_index["pop_size"]]
        @test length(brkga_data.previous[i].fitness) == param_values[param_index["pop_size"]]

        @test brkga_data.current[i].chromosomes == brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].chromosomes !== brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].fitness == brkga_data.previous[i].fitness
        @test brkga_data.current[i].fitness !== brkga_data.previous[i].fitness

        correct_order = true
        for j = 2:length(brkga_data.current[i].fitness)
            correct_order &= brkga_data.current[i].fitness[j-1] >=
                             brkga_data.current[i].fitness[j]
        end
        @test correct_order
    end

    @test brkga_data.initialized == true
    @test brkga_data.reset_phase == false

    param_values = copy(default_param_values)
    param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    for i = 1:brkga_data.num_independent_populations
        correct_order = true
        for j = 2:length(brkga_data.current[i].fitness)
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
    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)
    set_initial_population!(brkga_data, chromosomes)

    initialize!(brkga_data)

    for i = 1:brkga_data.num_independent_populations
        @test length(brkga_data.current) == param_values[param_index["num_independent_populations"]]
        @test length(brkga_data.current[i].chromosomes) == param_values[param_index["pop_size"]]
        @test length(brkga_data.current[i].fitness) == param_values[param_index["pop_size"]]

        @test length(brkga_data.previous) == param_values[param_index["num_independent_populations"]]
        @test length(brkga_data.previous[i].chromosomes) == param_values[param_index["pop_size"]]
        @test length(brkga_data.previous[i].fitness) == param_values[param_index["pop_size"]]

        @test brkga_data.current[i].chromosomes == brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].chromosomes !== brkga_data.previous[i].chromosomes
        @test brkga_data.current[i].fitness == brkga_data.previous[i].fitness
        @test brkga_data.current[i].fitness !== brkga_data.previous[i].fitness
    end

    decode!(chromosomes[1], instance)
    @test brkga_data.current[1].chromosomes[1] == chromosomes[1]

    # Create a local chromosome and applied the decoder on it.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    rand(local_rng, 1000)
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    decode!(local_chr, instance)

    # 4th chromosome must be the 1st generated due to the warmstart.
    @test brkga_data.current[1].chromosomes[4] == local_chr

    ########################
    # Test reset phase
    ########################

    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    # Create a local RNG and advance it until the same state as the internal
    # BrkgaData RNG after initialization.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    skip = brkga_data.num_independent_populations *
           brkga_data.population_size * brkga_data.chromosome_size
    rand(local_rng, 1000 + skip)

    # Assert the both generators are in the same state.
    @assert rand(brkga_data.rng) == rand(local_rng)

    brkga_data.reset_phase = true
    initialize!(brkga_data)

    # Create a local chromosome and applied the decoder on it.
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    decode!(local_chr, instance)

    @test brkga_data.current[1].chromosomes[1] == local_chr
end

################################################################################

@testset "reset!()" begin
    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)

    @test_throws ErrorException reset!(brkga_data)

    initialize!(brkga_data)

    # Create a local RNG and advance it until the same state as the internal
    # BrkgaData RNG after initialization.
    local_rng = MersenneTwister(param_values[param_index["seed"]])
    skip = brkga_data.num_independent_populations *
           brkga_data.population_size * brkga_data.chromosome_size

    # Assert the both generators are in the same state.
    rand(local_rng, 1000 + skip)
    @assert rand(brkga_data.rng) == rand(local_rng)

    # Create a local chromosome and applied the decoder on it.
    local_chr = rand(local_rng, param_values[param_index["chr_size"]])
    decode!(local_chr, instance)

    # Reset and test the first individual.
    reset!(brkga_data)
    @test brkga_data.current[1].chromosomes[1] == local_chr

    # After reset, the reset phase flag should be deactivated.
    @test brkga_data.reset_phase == false;
end
