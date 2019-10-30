################################################################################
# path_relink_tests.jl: unit tests for path relink routines of BrkgaMpIpr.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 06, 2018 by ceandrade
# Last update: Oct 30, 2019 by ceandrade
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

using Printf

include("util.jl")

next_pair(x::Int64) = (x + 1, x + 2)

function dist(x, y)
    total = 0.0
    for i in 1:length(x)
        if x[i] != y[i]
            total += 1.0
        end
    end
    return total
end

################################################################################

@testset "find_block_range()" begin
    SIZE::Int64 = 10
    for i in 1:SIZE
        @test BrkgaMpIpr.find_block_range(i, 1, SIZE) == i:i
    end

    @test BrkgaMpIpr.find_block_range(1, 2, SIZE) == 1:2
    @test BrkgaMpIpr.find_block_range(2, 2, SIZE) == 3:4
    @test BrkgaMpIpr.find_block_range(3, 2, SIZE) == 5:6
    @test BrkgaMpIpr.find_block_range(4, 2, SIZE) == 7:8
    @test BrkgaMpIpr.find_block_range(5, 2, SIZE) == 9:10

    @test BrkgaMpIpr.find_block_range(1, 3, SIZE) == 1:3
    @test BrkgaMpIpr.find_block_range(2, 3, SIZE) == 4:6
    @test BrkgaMpIpr.find_block_range(3, 3, SIZE) == 7:9
    @test BrkgaMpIpr.find_block_range(4, 3, SIZE) == 10:10
end

################################################################################

@testset "swap!()" begin
    a = [1, 2, 3]

    BrkgaMpIpr.swap!(a, 1, 3)
    @test a == [3, 2, 1]

    BrkgaMpIpr.swap!(a, 2, 3)
    @test a == [3, 1, 2]

    BrkgaMpIpr.swap!(a, 3, 1)
    @test a == [2, 1, 3]

    BrkgaMpIpr.swap!(a, 2, 1)
    @test a == [1, 2, 3]

    a = [1.1, 2.2, 3.3]
    BrkgaMpIpr.swap!(a, 2, 1)
    @test a == [2.2, 1.1, 3.3]

    a = ["a", "b", "c"]
    BrkgaMpIpr.swap!(a, 2, 1)
    @test a == ["b", "a", "c"]
end

################################################################################

@testset "|(::PathRelinkingResult, ::PathRelinkingResult)" begin
    @test (TOO_HOMOGENEOUS | TOO_HOMOGENEOUS) == TOO_HOMOGENEOUS
    @test (TOO_HOMOGENEOUS | NO_IMPROVEMENT) == NO_IMPROVEMENT
    @test (TOO_HOMOGENEOUS | ELITE_IMPROVEMENT) == ELITE_IMPROVEMENT
    @test (TOO_HOMOGENEOUS | BEST_IMPROVEMENT) == BEST_IMPROVEMENT
    @test (NO_IMPROVEMENT | NO_IMPROVEMENT) == NO_IMPROVEMENT
    @test (NO_IMPROVEMENT | ELITE_IMPROVEMENT) == ELITE_IMPROVEMENT
    @test (NO_IMPROVEMENT | BEST_IMPROVEMENT) == BEST_IMPROVEMENT
    @test (ELITE_IMPROVEMENT | ELITE_IMPROVEMENT) == ELITE_IMPROVEMENT
    @test (ELITE_IMPROVEMENT | BEST_IMPROVEMENT) == BEST_IMPROVEMENT
    @test (BEST_IMPROVEMENT | BEST_IMPROVEMENT) == BEST_IMPROVEMENT
end

################################################################################

@testset "direct_path_relink!()" begin
    start_time = time()

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MINIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = sum_decode!
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    ########################
    # Test fake homogeneity
    ########################
    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1234 #percentage::Float64
    )
    @test tmp == [Inf, Array{Float64, 1}()]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    param_values[param_index["opt_sense"]] = MAXIMIZE
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1234 #percentage::Float64
    )
    @test tmp == [-Inf, Array{Float64, 1}()]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Test maximum time
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["chr_size"]] = 10000
    param_values[param_index["instance"]] = Instance(10000)
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    start_time = time()
    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                5.0, #max_time::Float64,
                1.0 #percentage::Float64
    )

    # Test 5s. We must put a slack because it can fail on slow systems.
    @test ceil(time() - start_time) <= 10.0
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    start_time = time()
    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                10.0, #max_time::Float64,
                1.0 #percentage::Float64
    )

    # Test 10s. We must put a slack because it can fail on slow systems.
    @test ceil(time() - start_time) <= 15.0
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Test the relink
    ########################
    # **NOTE:** these tests may fail with the random number generation changes.
    # In such case, we have to figure out how to make this test better.

    data_path = joinpath(@__DIR__, "brkga_data_files")

    load_brkga_data(joinpath(data_path, "data_path_relink.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solutions_pr_direct.jld"))
    brkga_data.decode! = sum_decode!
    chr1 = 0
    chr2 = 0

    ###############
    # Block sizes
    ###############

    # Size 1
    chr1, chr2 = next_pair(chr2)
    block1 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block1[1] ≈ results["block1"][1]
    @test block1[2] ≈ results["block1"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 10
    chr1, chr2 = next_pair(chr2)
    block10 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block10[1] ≈ results["block10"][1]
    @test block10[2] ≈ results["block10"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 100
    chr1, chr2 = next_pair(chr2)
    block100 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                100, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block100[1] ≈ results["block100"][1]
    @test block100[2] ≈ results["block100"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 400
    chr1, chr2 = next_pair(chr2)
    block400 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                400, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block400[1] ≈ results["block400"][1]
    @test block400[2] ≈ results["block400"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 372
    chr1, chr2 = next_pair(chr2)
    block372 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                372, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block372[1] ≈ results["block372"][1]
    @test block372[2] ≈ results["block372"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ###############
    # Path sizes
    ###############

    # Path 10%
    chr1, chr2 = next_pair(chr2)
    path10 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1 #percentage::Float64
    )
    @test path10[1] ≈ results["path10"][1]
    @test path10[2] ≈ results["path10"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 30%
    chr1, chr2 = next_pair(chr2)
    path30 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.3 #percentage::Float64
    )
    @test path30[1] ≈ results["path30"][1]
    @test path30[2] ≈ results["path30"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 50%
    chr1, chr2 = next_pair(chr2)
    path50 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.50 #percentage::Float64
    )
    @test path50[1] ≈ results["path50"][1]
    @test path50[2] ≈ results["path50"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 100%
    chr1, chr2 = next_pair(chr2)
    path100 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                1.0 #percentage::Float64
    )
    @test path100[1] ≈ results["path100"][1]
    @test path100[2] ≈ results["path100"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ##############################
    # Simple distance function
    ##############################

    # x < y
    chr1, chr2 = next_pair(chr2)
    xy = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> x[1] < y[2], #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test xy[1] ≈ results["xy"][1]
    @test xy[2] ≈ results["xy"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # x > y
    chr1, chr2 = next_pair(chr2)
    yx = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> x[1] > y[2], #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test yx[1] ≈ results["yx"][1]
    @test yx[2] ≈ results["yx"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")
end

################################################################################

@testset "permutation_based_path_relink!()" begin
    start_time = time()

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = rank_decode!
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    ########################
    # Test fake homogeneity
    ########################

    brkga_data.current[1].chromosomes[1] = ones(param_values[param_index["chr_size"]])
    brkga_data.current[1].chromosomes[2] = ones(param_values[param_index["chr_size"]])

    tmp = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1234 #percentage::Float64
    )
    @test tmp == [-Inf, Array{Float64, 1}()]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    param_values[param_index["opt_sense"]] = MINIMIZE
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)
    brkga_data.current[1].chromosomes[1] = ones(param_values[param_index["chr_size"]])
    brkga_data.current[1].chromosomes[2] = ones(param_values[param_index["chr_size"]])

    tmp = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1234 #percentage::Float64
    )
    @test tmp == [Inf, Array{Float64, 1}()]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Test maximum time
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["chr_size"]] = 10000
    param_values[param_index["instance"]] = Instance(10000)
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    start_time = time()
    tmp = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                5.0, #max_time::Float64,
                1.0 #percentage::Float64
    )

    # Test 5s. We must put a slack because it can fail on slow systems.
    @test ceil(time() - start_time) <= 10.0
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    start_time = time()
    tmp = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                10.0, #max_time::Float64,
                1.0 #percentage::Float64
    )

    # Test 5s. We must put a slack because it can fail on slow systems.
    @test ceil(time() - start_time) <= 15.0
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    #########################
    # Test the relink
    #########################
    # **NOTE:** this test may fail with the random number generation changes.
    # In such case, we have to figure out how to make this test better.

    data_path = joinpath(@__DIR__, "brkga_data_files")

    load_brkga_data(joinpath(data_path, "data_path_relink.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solutions_pr_permutation_based.jld"))
    brkga_data.decode! = rank_decode!
    chr1 = 0
    chr2 = 0

    ###############
    # Block sizes
    ###############

    # Size 1
    chr1, chr2 = next_pair(chr2)
    block1 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 10
    chr1, chr2 = next_pair(chr2)
    block10 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block10[1] ≈ results["block10"][1]
    @test block10[2] ≈ results["block10"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 100
    chr1, chr2 = next_pair(chr2)
    block100 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                100, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block100[1] ≈ results["block100"][1]
    @test block100[2] ≈ results["block100"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 400
    chr1, chr2 = next_pair(chr2)
    block400 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                400, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block400[1] ≈ results["block400"][1]
    @test block400[2] ≈ results["block400"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Size 372
    chr1, chr2 = next_pair(chr2)
    block372 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                372, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test block372[1] ≈ results["block372"][1]
    @test block372[2] ≈ results["block372"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ###############
    # Path sizes
    ###############

    # Path 10%
    chr1, chr2 = next_pair(chr2)
    path10 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.1 #percentage::Float64
    )
    @test path10[1] ≈ results["path10"][1]
    @test path10[2] ≈ results["path10"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 30%
    chr1, chr2 = next_pair(chr2)
    path30 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.3 #percentage::Float64
    )
    @test path30[1] ≈ results["path30"][1]
    @test path30[2] ≈ results["path30"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 50%
    chr1, chr2 = next_pair(chr2)
    path50 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                0.5 #percentage::Float64
    )
    @test path50[1] ≈ results["path50"][1]
    @test path50[2] ≈ results["path50"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # Path 100%
    chr1, chr2 = next_pair(chr2)
    path100 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120.0, #max_time::Float64,
                1.0 #percentage::Float64
    )
    @test path100[1] ≈ results["path100"][1]
    @test path100[2] ≈ results["path100"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ##############################
    # Simple distance function
    ##############################

    # x < y
    chr1, chr2 = next_pair(chr2)
    xy = BrkgaMpIpr.permutation_based_path_relink!(
            brkga_data, #brkga_data::BrkgaData,
            brkga_data.current[1].chromosomes[chr1], #chromosome1
            brkga_data.current[1].chromosomes[chr2], #chromosome2
            (x, y) -> x[1] < y[2], #affect_solution::Function,
            10, #block_size::Int64,
            120.0, #max_time::Float64,
            0.5 #percentage::Float64
    )
    @test xy[1] ≈ results["xy"][1]
    @test xy[2] ≈ results["xy"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    # x > y
    chr1, chr2 = next_pair(chr2)
    yx = BrkgaMpIpr.permutation_based_path_relink!(
            brkga_data, #brkga_data::BrkgaData,
            brkga_data.current[1].chromosomes[chr1], #chromosome1
            brkga_data.current[1].chromosomes[chr2], #chromosome2
            (x, y) -> x[1] > y[2], #affect_solution::Function,
            10, #block_size::Int64,
            120.0, #max_time::Float64,
            0.5 #percentage::Float64
    )
    @test yx[1] ≈ results["yx"][1]
    @test yx[2] ≈ results["yx"][2]
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")
end

################################################################################

@testset "path_relink!()" begin
    start_time = time()

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_independent_populations = 1
    brkga_data = build_brkga(param_values...)

    ########################
    # Test arguments
    ########################

    @test_throws ErrorException path_relink!(brkga_data, #::BrkgaData,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            1, #block_size::Int64,
                                            120.0, #max_time::Float64,
                                            1.0 #percentage::Float64
                                            )

    initialize!(brkga_data)

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            1, #block_size::Int64,
                                            120.0, #max_time::Float64,
                                            -0.1 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            1, #block_size::Int64,
                                            120.0, #max_time::Float64,
                                            1.01 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            0, #block_size::Int64,
                                            120.0, #max_time::Float64,
                                            1.0 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            -10, #block_size::Int64,
                                            120.0, #max_time::Float64,
                                            1.0 #percentage::Float64
                                            )

    ########################
    # Test fake homogeneity
    ########################

    @test TOO_HOMOGENEOUS == path_relink!(
        brkga_data, #::BrkgaData,
        DIRECT, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        (x, y) -> 0.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        10, #number_pairs::Int64,
        1.0, #minimum_distance::Float64,
        1, #block_size::Int64,
        2.0, #max_time::Float64,
        1.0 #percentage::Float64
    )
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    #######################
    # Test the direct path relink
    #######################

    # Let's save the data to use after this.
    original_brkga_data = deepcopy(brkga_data)

    @test BEST_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        DIRECT, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        (x, y) -> 1.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        10, #number_pairs::Int64,
        0.5, #minimum_distance::Float64,
        1, #block_size::Int64,
        10.0, #max_time::Float64,
        1.0, #percentage::Float64
    )

    # Do not change the best.
    @test get_best_fitness(brkga_data) > get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) != get_best_chromosome(original_brkga_data)
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Test the permutation path relink
    ########################

    ########################
    # No improvement found
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = sum_decode!
    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_independent_populations = 1
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    # Create a homogeneous population and instance, so that any combination
    # doesn't generate improvements.
    brkga_data.problem_instance.data .=
        ones(length(brkga_data.problem_instance.data))

    population = brkga_data.current[1]
    for chr in population.chromosomes
        chr .= ones(length(chr))
    end

    # Let's re-decode everything.
    for i in 1:brkga_data.params.population_size
        value = brkga_data.decode!(population.chromosomes[i],
                                   brkga_data.problem_instance, false)
        population.fitness[i] = (value, i)
    end

    original_brkga_data = deepcopy(brkga_data)

    @test NO_IMPROVEMENT == path_relink!(
            brkga_data, #::BrkgaData,
            DIRECT, #::PathRelinkingType,
            BESTSOLUTION, # PathRelinkingSelection
            (x, y) -> 1.0, #compute_distance::Function,
            (x, y) -> true, #affect_solution::Function,
            10, #number_pairs::Int64,
            0.5, #minimum_distance::Float64,
            1, #block_size::Int64,
            10.0, #max_time::Float64,
            1.0, #percentage::Float64
        )

    # Do not change at all.
    @test get_best_fitness(brkga_data) == get_best_fitness(original_brkga_data)
    @test get_current_population(brkga_data, 1) !=
          get_current_population(original_brkga_data, 1)
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Improvement found, but not in the best solution.
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 100
    param_values[param_index["instance"]] = Instance(100)
    param_values[param_index["decode!"]] = rank_decode!
    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_independent_populations = 1
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)
    original_brkga_data = deepcopy(brkga_data)

    @test ELITE_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        PERMUTATION, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        (x, y) -> 1.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        0, #number_pairs::Int64,
        0.5, #minimum_distance::Float64,
        1, #block_size::Int64,
        0.0, #max_time::Float64,
        1.0, #percentage::Float64
    )

    # Do not change the best.
    @test get_best_fitness(brkga_data) ≈ get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) ≈ get_best_chromosome(original_brkga_data)

    # But population changes.
    @test get_current_population(brkga_data, 1) !=
          get_current_population(original_brkga_data, 1)
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Improvement found
    ########################
    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = rank_decode!
    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_independent_populations = 1
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    population = brkga_data.current[1]

    # First, create a homogeneous population and instance
    brkga_data.problem_instance.data .=
        zeros(length(brkga_data.problem_instance.data))
    for chr in population.chromosomes
        chr .= zeros(length(chr))
    end

    # # Now, inject two chromossomes with partial permutations. These chromosomes
    # should be selected during the path relink.
    half = div(param_values[param_index["chr_size"]], 2)
    for i in 1:half
        population.chromosomes[1][i] = i / 10.0
        population.chromosomes[2][i + half] = (i + half) / 10.0
    end

    # Let's re-decode everything.
    for i in 1:brkga_data.params.population_size
        value = brkga_data.decode!(population.chromosomes[i],
                                   brkga_data.problem_instance, false)
        population.fitness[i] = (value, i)
    end

    original_brkga_data = deepcopy(brkga_data)

    @test BEST_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        PERMUTATION, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        dist, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        0, #number_pairs::Int64,
        1.0, #minimum_distance::Float64,
        1, #block_size::Int64,
        0.0, #max_time::Float64,
        1.0, #percentage::Float64
    )

    @test get_best_fitness(brkga_data) > get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) != get_best_chromosome(original_brkga_data)
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")

    ########################
    # Five populations
    ########################

    param_values = deepcopy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MINIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = sum_decode!
    brkga_params = param_values[param_index["brkga_params"]]
    brkga_params.num_independent_populations = 5
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    # First, create a homogeneous population and instance
    brkga_data.problem_instance.data .=
        ones(length(brkga_data.problem_instance.data))

    for pop in brkga_data.current
        for chr in pop.chromosomes
            chr .= zeros(length(chr))
        end
    end

    # Keep populations 1, 2, and 5. Change 3 and 4.
    brkga_data.current[3].chromosomes[1] =
        rand(brkga_data.rng, brkga_data.chromosome_size)
    brkga_data.current[4].chromosomes[1] =
        rand(brkga_data.rng, brkga_data.chromosome_size)

    # Let's re-decode everything.
    for pop in brkga_data.current
        for i in 1:brkga_data.params.population_size
            value = brkga_data.decode!(pop.chromosomes[i],
                                       brkga_data.problem_instance, false)
            pop.fitness[i] = (value, i)
        end
    end

    original_brkga_data = deepcopy(brkga_data)

    @test BEST_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        DIRECT, #::PathRelinkingType,
        BESTSOLUTION, # PathRelinkingSelection
        dist, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        0, #number_pairs::Int64,
        1.0, #minimum_distance::Float64,
        1, #block_size::Int64,
        0.0, #max_time::Float64,
        1.0, #percentage::Float64
    )

    @test get_best_fitness(brkga_data) < get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) != get_best_chromosome(original_brkga_data)
    println("Elapsed time: $(@sprintf("%.2f", time() - start_time))")
end
