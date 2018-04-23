################################################################################
# evolution_tests.jl: unit tests for evolutionary routines of BrkgaMpIpr.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Apr 19, 2018 by ceandrade
# Last update: Apr 23, 2018 by ceandrade
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

function load_brkga_data(filename::String, brkga_data::BrkgaData)
    tmp = load(filename)

    brkga_data.opt_sense = tmp["opt_sense"]
    brkga_data.chromosome_size = tmp["chromosome_size"]
    brkga_data.population_size = tmp["population_size"]
    brkga_data.elite_size = tmp["elite_size"]
    brkga_data.num_mutants = tmp["num_mutants"]
    brkga_data.num_elite_parents = tmp["num_elite_parents"]
    brkga_data.total_parents = tmp["total_parents"]
    brkga_data.bias = tmp["bias"]
    brkga_data.num_independent_populations = tmp["num_independent_populations"]
    brkga_data.evolutionary_mechanism_on = tmp["evolutionary_mechanism_on"]

    # TODO (ceandrade): list the path relink parameters here.

    # FIXME (ceandrade): the following doesn't work because it tries to
    # load the decoder function from the file. So, we rebuild the instance.
    # brkga_data.problem_instance = tmp["problem_instance"],
    brkga_data.problem_instance =
        TestInstance.Instance(brkga_data.chromosome_size)

    # NOTE (ceandrade): currently, JLD cannot save functions.
    brkga_data.decode! = decode!

    brkga_data.rng = tmp["rng"]
    brkga_data.previous = tmp["previous"]
    brkga_data.current = tmp["current"]
    brkga_data.total_bias_weight = tmp["total_bias_weight"]
    brkga_data.shuffled_individuals = tmp["shuffled_individuals"]
    brkga_data.parents_ordered = tmp["parents_ordered"]
    brkga_data.initialized = tmp["initialized"]
    brkga_data.reset_phase = tmp["reset_phase"]

    # NOTE (ceandrade): currently, JLD cannot save functions.
    if brkga_data.bias == LOGINVERSE
        set_bias_custom_function!(brkga_data, r -> 1.0 / log1p(r))
    elseif brkga_data.bias == LINEAR
        set_bias_custom_function!(brkga_data, r -> 1.0 / r)
    elseif brkga_data.bias == QUADRATIC
        set_bias_custom_function!(brkga_data, r -> r ^ -2.0)
    elseif brkga_data.bias == CUBIC
        set_bias_custom_function!(brkga_data, r -> r ^ -3.0)
    elseif brkga_data.bias == EXPONENTIAL
        set_bias_custom_function!(brkga_data, r -> exp(-r))
    elseif brkga_data.bias == CONSTANT
        set_bias_custom_function!(brkga_data, (::Int64) -> 1.0 / total_parents)
    end

    nothing
end

################################################################################

@testset "evolve_population!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["num_elite_parents"]] = 2
    param_values[param_index["total_parents"]] = 3
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

    for i in 1:brkga_data.num_independent_populations
        evolve_population!(brkga_data, i)

        @test current[i].chromosomes == brkga_data.previous[i].chromosomes
        @test current[i].fitness == brkga_data.previous[i].fitness

        @test previous[i].chromosomes != brkga_data.current[i].chromosomes
        @test previous[i].fitness != brkga_data.current[i].fitness
    end

    ########################
    # Test the evolutionary mechanism
    ########################
    # **NOTE:** this test may fail with the random number generation changes.
    # In such case, we have to figure out how to make this test better.

    data_path = joinpath(@__DIR__, "brkga_data_files")

    ########################
    # Data 1
    load_brkga_data(joinpath(data_path, "data1.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solution1.jld"))

    evolve_population!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness1"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome1"]

    evolve_population!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness2"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome2"]

    for _ in 1:100
        evolve_population!(brkga_data, 1)
    end
    @test get_best_fitness(brkga_data) ≈ results["fitness102"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome102"]

    ########################
    # Data 2
    load_brkga_data(joinpath(data_path, "data2.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solution2.jld"))

    evolve_population!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness1"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome1"]

    evolve_population!(brkga_data, 2)
    @test get_best_fitness(brkga_data) ≈ results["fitness2"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome2"]

    for _ in 1:100
        evolve_population!(brkga_data, 1)
        evolve_population!(brkga_data, 2)
    end
    @test get_best_fitness(brkga_data) ≈ results["fitness102"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome102"]

    ########################
    # Data 3
    load_brkga_data(joinpath(data_path, "data3.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solution3.jld"))

    evolve_population!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness1"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome1"]

    for _ in 2:brkga_data.num_independent_populations
        evolve_population!(brkga_data, 2)
    end
    @test get_best_fitness(brkga_data) ≈ results["fitness2"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome2"]

    for _ in 1:100
        for i in 1:brkga_data.num_independent_populations
            evolve_population!(brkga_data, i)
        end
    end
    @test get_best_fitness(brkga_data) ≈ results["fitness102"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome102"]

    ########################
    # Data 4 (traditional BRKGA)
    load_brkga_data(joinpath(data_path, "data4.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solution4.jld"))

    rho = 0.75
    set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)

    evolve_population!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness1"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome1"]

    evolve_population!(brkga_data, 2)
    evolve_population!(brkga_data, 3)
    @test get_best_fitness(brkga_data) ≈ results["fitness2"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome2"]

    for _ in 1:100
        for i in 1:brkga_data.num_independent_populations
            evolve_population!(brkga_data, i)
        end
    end
    @test get_best_fitness(brkga_data) ≈ results["fitness102"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome102"]
end

################################################################################

@testset "evolve!()" begin
    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)

    # Not initialized
    @test_throws ErrorException evolve!(brkga_data)

    initialize!(brkga_data)

    @test_throws ArgumentError evolve!(brkga_data, 0)
    @test_throws ArgumentError evolve!(brkga_data, -10)

    #################################
    # Several evolutionary iterations
    #################################

    data_path = joinpath(@__DIR__, "brkga_data_files")
    load_brkga_data(joinpath(data_path, "data5.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solution5.jld"))

    evolve!(brkga_data, 1)
    @test get_best_fitness(brkga_data) ≈ results["fitness1"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome1"]

    evolve!(brkga_data, 10)
    @test get_best_fitness(brkga_data) ≈ results["fitness2"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome2"]

    evolve!(brkga_data, 100)
    @test get_best_fitness(brkga_data) ≈ results["fitness102"]
    @test get_best_chromosome(brkga_data) ≈ results["chromosome102"]
end
