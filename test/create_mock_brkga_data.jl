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
# Last update: Jun 11, 2018 by ceandrade
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

include("TestInstance.jl")
include("TestDecoder.jl")

using BrkgaMpIpr
using JLD
using TestDecoder
using TestInstance

# Makes easy to change specific position on the parameters vector below.
const param_names = ["instance", "decode!", "opt_sense", "seed", "chr_size",
                     "pop_size", "elite_percentage", "mutants_percentage",
                     "evolutionary_mechanism_on", "num_elite_parents",
                     "total_parents", "bias", "num_independent_populations"]

# Reverse index.
const param_index = Dict([v => i for (i, v) in enumerate(param_names)])

# Holds the parameters to build new BrkgaData.
param_values = Array{Any, 1}(length(param_names))

################################################################################
# Write file
################################################################################

function write_data(filename::String, data::BrkgaData)
    save(filename,
        "opt_sense", brkga_data.opt_sense,
        "chromosome_size", brkga_data.chromosome_size,
        "population_size", brkga_data.population_size,
        "elite_size", brkga_data.elite_size,
        "num_mutants", brkga_data.num_mutants,
        "num_elite_parents", brkga_data.num_elite_parents,
        "total_parents", brkga_data.total_parents,
        "bias", brkga_data.bias,
        "num_independent_populations", brkga_data.num_independent_populations,
        "evolutionary_mechanism_on", brkga_data.evolutionary_mechanism_on,
        # TODO (ceandrade): list the path relink parameters here.
        "problem_instance", brkga_data.problem_instance,
        # NOTE (ceandrade): currently, JLD cannot save functions.
        # decode!::Function
        "rng", brkga_data.rng,
        "previous", brkga_data.previous,
        "current", brkga_data.current,
        # NOTE (ceandrade): currently, JLD cannot save functions.
        # "bias_function", brkga_data.bias_function,
        "total_bias_weight", brkga_data.total_bias_weight,
        "shuffled_individuals", brkga_data.shuffled_individuals,
        "parents_ordered", brkga_data.parents_ordered,
        "initialized", brkga_data.initialized,
        "reset_phase", brkga_data.reset_phase
    )
end

################################################################################
# Configuration 1
################################################################################

chromosome_size = 100
instance = Instance(chromosome_size)
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = decode!
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

print("\n\n> Building configuration 1")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 1")
write_data("brkga_data_files/data1.jld", brkga_data)

print("\n> Evolving population 1")

evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

evolve_population!(brkga_data, 1)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    evolve_population!(brkga_data, 1)
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save("brkga_data_files/best_solution1.jld",
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
param_values[param_index["decode!"]] = decode!
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

print("\n\n> Building configuration 2")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 2")
write_data("brkga_data_files/data2.jld", brkga_data)

print("\n> Evolving population 2")

evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

evolve_population!(brkga_data, 2)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    evolve_population!(brkga_data, 1)
    evolve_population!(brkga_data, 2)
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save("brkga_data_files/best_solution2.jld",
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
param_values[param_index["decode!"]] = decode!
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

print("\n\n> Building configuration 3")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing configuration 3")
write_data("brkga_data_files/data3.jld", brkga_data)

print("\n> Evolving population 3")

evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

for _ in 2:brkga_data.num_independent_populations
    evolve_population!(brkga_data, 2)
end

fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    for i in 1:brkga_data.num_independent_populations
        evolve_population!(brkga_data, i)
    end
end
fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save("brkga_data_files/best_solution3.jld",
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
param_values[param_index["decode!"]] = decode!
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

print("\n\n> Building configuration 4")
brkga_data = build_brkga(param_values...)

rho = 0.75
set_bias_custom_function!(brkga_data, x -> x ≈ 1.0 ? rho : 1.0 - rho)
initialize!(brkga_data)

print("\n> Writing configuration 4")
write_data("brkga_data_files/data4.jld", brkga_data)

print("\n> Evolving population 4")

evolve_population!(brkga_data, 1)
fitness1 = get_best_fitness(brkga_data)
chromosome1 = get_best_chromosome(brkga_data)

evolve_population!(brkga_data, 2)
evolve_population!(brkga_data, 3)
fitness2 = get_best_fitness(brkga_data)
chromosome2 = get_best_chromosome(brkga_data)

for _ in 1:100
    for i in 1:brkga_data.num_independent_populations
        evolve_population!(brkga_data, i)
    end
end

fitness102 = get_best_fitness(brkga_data)
chromosome102 = get_best_chromosome(brkga_data)

save("brkga_data_files/best_solution4.jld",
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
param_values[param_index["decode!"]] = decode!
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

save("brkga_data_files/best_solution5.jld",
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
param_values[param_index["decode!"]] = decode!
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

print("\n\n> Building configuration for path relink")
brkga_data = build_brkga(param_values...)
initialize!(brkga_data)

print("\n> Writing path relink")
write_data("brkga_data_files/data_path_relink.jld", brkga_data)

# Create some test data for path relink methods.
print("\n> Path relinking population ")

###############
# Block sizes
###############

# Size 1
block1 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           1, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

# Size 10
block10 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           1, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

# Size 100
block100 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           100, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

# Size 400
block400 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           400, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

# Size 372
block372 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           372, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

###############
# Path sizes
###############

# Path 10%
path10 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           0.1 #percentage::Float64
                           )

# Path 30%
path30 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           0.3 #percentage::Float64
                           )

# Path 50%
path50 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           0.50 #percentage::Float64
                           )

# Path 100%
path100 = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> true, #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           1.0 #percentage::Float64
                           )

##############################
# Simple distance function
##############################

# x < y
xy = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> x[1] < y[2], #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

# x > y
yx = direct_path_relink!(brkga_data, #brkga_data::BrkgaData,
                           1, #population_index::Int64,
                           1, #chr1_index::Int64,
                           2, #chr2_index::Int64,
                           (x, y) -> x[1] > y[2], #distance_function::Function,
                           10, #block_size::Int64,
                           120, #max_time::Int64,
                           0.5 #percentage::Float64
                           )

save("brkga_data_files/best_solutions_direct_pr.jld",
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