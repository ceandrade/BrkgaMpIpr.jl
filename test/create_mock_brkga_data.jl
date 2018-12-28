#!/usr/bin/env julia
################################################################################
# create_mock_brkga_data.jl: create BRKGA data mock object to be used in
# the evolutionary tests.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Apr 20, 2018 by ceandrade
# Last update: Dec 28, 2018 by ceandrade
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

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using BrkgaMpIpr

include("test_instance.jl")
include("test_decoders.jl")
include("util.jl")

# Makes easy to change specific position on the parameters vector below.
const param_names = ["instance", "decode!", "opt_sense", "seed", "chr_size",
                     "pop_size", "elite_percentage", "mutants_percentage",
                     "evolutionary_mechanism_on", "num_elite_parents",
                     "total_parents", "bias", "num_independent_populations",
                     "pr_number_pairs", "pr_minimum_distance", "pr_type",
                     "pr_selection", "alpha_block_size", "pr_percentage"]

# Reverse index.
const param_index = Dict([v => i for (i, v) in enumerate(param_names)])

# Holds the parameters to build new BrkgaData.
param_values = Array{Any, 1}(undef, length(param_names))

################################################################################
# Configuration 1
################################################################################

chromosome_size = 100
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MAXIMIZE
param_values[param_index["seed"]] = 3979164113692134205
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 100
param_values[param_index["elite_percentage"]] = 0.3
param_values[param_index["mutants_percentage"]] = 0.1
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 1
param_values[param_index["total_parents"]] = 2
param_values[param_index["bias"]] = LOGINVERSE
param_values[param_index["num_independent_populations"]] = 1
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration 1")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 1")
write_data("brkga_data_files/data1.jld", brkga_data)

print("\n> Evolving population 1")

BrkgaMpIpr.evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

BrkgaMpIpr.evolve_population!(brkga_data, 1)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    BrkgaMpIpr.evolve_population!(brkga_data, 1)
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save(File(format"JLD", "brkga_data_files/best_solution1.jld"),
    "fitness1", fitness1,
    "chromosome1", chromosome1,
    "fitness2", fitness2,
    "chromosome2", chromosome2,
    "fitness102", fitness102,
    "chromosome102", chromosome102,
)

################################################################################
# Configuration 2
################################################################################

chromosome_size = 1000
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MINIMIZE
param_values[param_index["seed"]] = 1297832326904308
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 500
param_values[param_index["elite_percentage"]] = 0.25
param_values[param_index["mutants_percentage"]] = 0.25
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 5
param_values[param_index["total_parents"]] = 50
param_values[param_index["bias"]] = QUADRATIC
param_values[param_index["num_independent_populations"]] = 2
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration 2")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 2")
write_data("brkga_data_files/data2.jld", brkga_data)

print("\n> Evolving population 2")

BrkgaMpIpr.evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

BrkgaMpIpr.evolve_population!(brkga_data, 2)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    BrkgaMpIpr.evolve_population!(brkga_data, 1)
    BrkgaMpIpr.evolve_population!(brkga_data, 2)
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save(File(format"JLD", "brkga_data_files/best_solution2.jld"),
    "fitness1", fitness1,
    "chromosome1", chromosome1,
    "fitness2", fitness2,
    "chromosome2", chromosome2,
    "fitness102", fitness102,
    "chromosome102", chromosome102,
)

################################################################################
# Configuration 3
################################################################################

chromosome_size = 500
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MINIMIZE
param_values[param_index["seed"]] = 2536246074066680359
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 100
param_values[param_index["elite_percentage"]] = 0.35
param_values[param_index["mutants_percentage"]] = 0.17
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 3
param_values[param_index["total_parents"]] = 5
param_values[param_index["bias"]] = EXPONENTIAL
param_values[param_index["num_independent_populations"]] = 5
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration 3")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 3")
write_data("brkga_data_files/data3.jld", brkga_data)

print("\n> Evolving population 3")

BrkgaMpIpr.evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

for _ in 2:brkga_data.num_independent_populations
    BrkgaMpIpr.evolve_population!(brkga_data, 2)
end

fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    for i in 1:brkga_data.num_independent_populations
        BrkgaMpIpr.evolve_population!(brkga_data, i)
    end
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save(File(format"JLD", "brkga_data_files/best_solution3.jld"),
    "fitness1", fitness1,
    "chromosome1", chromosome1,
    "fitness2", fitness2,
    "chromosome2", chromosome2,
    "fitness102", fitness102,
    "chromosome102", chromosome102,
)

################################################################################
# Configuration 4 (traditional BRKGA)
################################################################################

chromosome_size = 500
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MINIMIZE
param_values[param_index["seed"]] = 2947804214766190222
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 100
param_values[param_index["elite_percentage"]] = 0.30
param_values[param_index["mutants_percentage"]] = 0.15
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 1
param_values[param_index["total_parents"]] = 2
param_values[param_index["bias"]] = LOGINVERSE
param_values[param_index["num_independent_populations"]] = 3
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration 4")
brkga_data = build_brkga(param_values...)

rho = 0.75
set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)
initialize!(brkga_data)

print("\n> Writing configuration 4")
write_data("brkga_data_files/data4.jld", brkga_data)

print("\n> Evolving population 4")

BrkgaMpIpr.evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

BrkgaMpIpr.evolve_population!(brkga_data, 2)
BrkgaMpIpr.evolve_population!(brkga_data, 3)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    for i in 1:brkga_data.num_independent_populations
        BrkgaMpIpr.evolve_population!(brkga_data, i)
    end
end

fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save(File(format"JLD", "brkga_data_files/best_solution4.jld"),
    "fitness1", fitness1,
    "chromosome1", chromosome1,
    "fitness2", fitness2,
    "chromosome2", chromosome2,
    "fitness102", fitness102,
    "chromosome102", chromosome102,
)

################################################################################
# Configuration 5 for evolve!
################################################################################

chromosome_size = 100
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MINIMIZE
param_values[param_index["seed"]] = 4659930950303615118
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 100
param_values[param_index["elite_percentage"]] = 0.30
param_values[param_index["mutants_percentage"]] = 0.20
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 2
param_values[param_index["total_parents"]] = 3
param_values[param_index["bias"]] = LOGINVERSE
param_values[param_index["num_independent_populations"]] = 3
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration 5")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 5")
write_data("brkga_data_files/data5.jld", brkga_data)

print("\n> Evolving population 5")

evolve!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

evolve!(brkga_data, 10)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

evolve!(brkga_data, 100)
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save(File(format"JLD", "brkga_data_files/best_solution5.jld"),
    "fitness1", fitness1,
    "chromosome1", chromosome1,
    "fitness2", fitness2,
    "chromosome2", chromosome2,
    "fitness102", fitness102,
    "chromosome102", chromosome102,
)

################################################################################
# Configuration for path relink
################################################################################

chromosome_size = 1000
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = sum_decode!
param_values[param_index["opt_sense"]] = MAXIMIZE
param_values[param_index["seed"]] = 16986526969459
param_values[param_index["chr_size"]] = chromosome_size
param_values[param_index["pop_size"]] = 100
param_values[param_index["elite_percentage"]] = 0.25
param_values[param_index["mutants_percentage"]] = 0.25
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 2
param_values[param_index["total_parents"]] = 3
param_values[param_index["bias"]] = QUADRATIC
param_values[param_index["num_independent_populations"]] = 1
param_values[param_index["pr_number_pairs"]] = 0
param_values[param_index["pr_minimum_distance"]] = 0.0
param_values[param_index["pr_type"]] = DIRECT
param_values[param_index["pr_selection"]] = BESTSOLUTION
param_values[param_index["alpha_block_size"]] = 1.0
param_values[param_index["pr_percentage"]] = 1.0

print("\n\n> Building configuration for path relink")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing path relink")
write_data("brkga_data_files/data_path_relink.jld", brkga_data)

next_pair(x::Int64) = (x + 1, x + 2)

# Create some test data for path relink methods.
for (func, decoder, name) in
    [(BrkgaMpIpr.direct_path_relink!, sum_decode!, "direct"),
     (BrkgaMpIpr.permutation_based_path_relink!, rank_decode!, "permutation_based")]

    print("\n> Generating results for tests for " * name)
    brkga_data.decode! = decoder

    chr1 = 0
    chr2 = 0

    ###############
    # Block sizes
    ###############

    # Size 1
    chr1, chr2 = next_pair(chr2)
    block1 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Size 10
    chr1, chr2 = next_pair(chr2)
    block10 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Size 100
    chr1, chr2 = next_pair(chr2)
    block100 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                100, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Size 400
    chr1, chr2 = next_pair(chr2)
    block400 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                400, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Size 372
    chr1, chr2 = next_pair(chr2)
    block372 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                372, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    ###############
    # Path sizes
    ###############

    # Path 10%
    chr1, chr2 = next_pair(chr2)
    path10 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.1 #percentage::Float64
    )

    # Path 30%
    chr1, chr2 = next_pair(chr2)
    path30 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.3 #percentage::Float64
    )

    # Path 50%
    chr1, chr2 = next_pair(chr2)
    path50 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Path 100%
    chr1, chr2 = next_pair(chr2)
    path100 = func(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                1.0 #percentage::Float64
    )

    ##############################
    # Simple distance function
    ##############################

    # x < y
    chr1, chr2 = next_pair(chr2)
    xy = func(
            brkga_data, #brkga_data::BrkgaData,
            brkga_data.current[1].chromosomes[chr1], #chromosome1
            brkga_data.current[1].chromosomes[chr2], #chromosome2
            (x, y) -> x[1] < y[2], #affect_solution::Function,
            10, #block_size::Int64,
            120, #max_time::Int64,
            0.5 #percentage::Float64
    )

    # x > y
    chr1, chr2 = next_pair(chr2)
    yx = func(
            brkga_data, #brkga_data::BrkgaData,
            brkga_data.current[1].chromosomes[chr1], #chromosome1
            brkga_data.current[1].chromosomes[chr2], #chromosome2
            (x, y) -> x[1] > y[2], #affect_solution::Function,
            10, #block_size::Int64,
            120, #max_time::Int64,
            0.5 #percentage::Float64
    )

    save(File(format"JLD", "brkga_data_files/best_solutions_pr_" *
                           name * ".jld"),
         "block1", block1,
         "block10", block10,
         "block100", block100,
         "block400", block400,
         "block372", block372,
         "path10", path10,
         "path30", path30,
         "path50", path50,
         "path100", path100,
         "xy", xy,
         "yx", yx
    )
end
