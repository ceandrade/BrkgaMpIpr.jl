################################################################################
# types_tests.jl: unit tests for type handling.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Oct, 30 2018 by ceandrade
# Last update: Oct, 30 2018 by ceandrade
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

@testset "Population constructors" begin
    population = BrkgaMpIpr.Population()
    @test population.chromosomes == Array{Array{Float64, 1}, 1}()
    @test population.fitness == Array{Tuple{Float64, Int64}, 1}()

    rng = Random.MersenneTwister(270000001)

    population_size = 100
    chromosome_size = 100
    population.fitness =
        Array{Tuple{Float64, Int64}, 1}(undef, population_size)

    for i in 1:population_size
        push!(population.chromosomes, rand(rng, chromosome_size))
        population.fitness[i] = rand(rng), i
    end

    population2 = BrkgaMpIpr.Population(population)

    @test population2.chromosomes == population.chromosomes
    @test population2.fitness == population.fitness

    population2.chromosomes[1] = rand(rng, chromosome_size)
    population2.fitness[1] = rand(rng), population_size + 1

    @test population2.chromosomes != population.chromosomes
    @test population2.fitness != population.fitness
end

################################################################################

@testset "BrkgaParams constructors" begin
    brkga_params = BrkgaMpIpr.BrkgaParams()

    @test brkga_params.population_size == 0
    @test brkga_params.elite_percentage == 0.0
    @test brkga_params.mutants_percentage == 0.0
    @test brkga_params.num_elite_parents == 0
    @test brkga_params.total_parents == 0
    @test brkga_params.bias_type == CONSTANT
    @test brkga_params.num_independent_populations == 0
    @test brkga_params.pr_number_pairs == 0
    @test brkga_params.pr_minimum_distance == 0.0
    @test brkga_params.pr_type == DIRECT
    @test brkga_params.pr_selection == BESTSOLUTION
    @test brkga_params.alpha_block_size == 0.0
    @test brkga_params.pr_percentage == 0.0
end

################################################################################

@testset "ExternalControlParams constructors" begin
    extra_params = BrkgaMpIpr.ExternalControlParams()

    @test extra_params.exchange_interval == 0
    @test extra_params.num_exchange_indivuduals == 0.0
    @test extra_params.reset_interval == 0.0
end

