################################################################################
# evolution_tests.jl: unit tests for evolutionary routines of BrkgaMpIpr.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Apr 19, 2018 by ceandrade
# Last update: Apr 20, 2018 by ceandrade
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

# @testset "evolve!()" begin
#     param_values = copy(default_param_values)
#     brkga_data = build_brkga(param_values...)
#
#     # Not initialized
#     @test_throws ErrorException evolve!(brkga_data)
#
#     initialize!(brkga_data)
#
#     @test_throws ArgumentError evolve!(brkga_data, 0)
#     @test_throws ArgumentError evolve!(brkga_data, -10)
#
#     evolve!(brkga_data, 2)
#
#     # # Create a local RNG and advance it until the same state as the internal
#     # # BrkgaData RNG after initialization.
#     # local_rng = MersenneTwister(param_values[param_index["seed"]])
#     # skip = brkga_data.num_independent_populations *
#     #        brkga_data.population_size * brkga_data.chromosome_size
#     #
#     # # Assert the both generators are in the same state.
#     # rand(local_rng, 1000 + skip)
#     # @assert rand(brkga_data.rng) == rand(local_rng)
#     #
#     # # Create a local chromosome and applied the decoder on it.
#     # local_chr = rand(local_rng, param_values[param_index["chr_size"]])
#     # decode!(local_chr, instance)
#     #
#     # # Reset and test the first individual.
#     # reset!(brkga_data)
#     # @test brkga_data.current[1].chromosomes[1] == local_chr
#     #
#     # # After reset, the reset phase flag should be deactivated.
#     # @test brkga_data.reset_phase == false;
# end

################################################################################

@testset "evolve_population!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["num_elite_parents"]] = 2
    param_values[param_index["total_parents"]] = 3
    # param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)

    # Not initialized
    @test_throws ErrorException evolve_population!(brkga_data, 1)

    initialize!(brkga_data)

    @test_throws ArgumentError evolve_population!(brkga_data, 0)
    @test_throws ArgumentError evolve_population!(brkga_data, -1)
    @test_throws ArgumentError evolve_population!(brkga_data,
        brkga_data.num_independent_populations + 1)

    # Save previous and current populations locally
    previous = deepcopy(brkga_data.previous)
    current = deepcopy(brkga_data.current)

    ########################
    # Test if algorithm swaps the populations correctly
    ########################

    evolve_population!(brkga_data, 1)

    @test current[1].chromosomes == brkga_data.previous[1].chromosomes
    @test current[1].fitness == brkga_data.previous[1].fitness

    @test previous[1].chromosomes != brkga_data.current[1].chromosomes
    @test previous[1].fitness != brkga_data.current[1].fitness

    ########################
    # Test if the best solution has change
    ########################
    # **NOTE:** this test may fail with the random number generation changes.
    # In such case, we have to figure out how to make this test better.
end
