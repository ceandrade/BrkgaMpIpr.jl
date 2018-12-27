################################################################################
# path_relink_tests.jl: unit tests for path relink routines of BrkgaMpIpr.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 06, 2018 by ceandrade
# Last update: Dec 26, 2018 by ceandrade
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

include("util.jl")

next_pair(x::Int64) = (x + 1, x + 2)

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

@testset "direct_path_relink!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MINIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["decode!"]] = decode!
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
                120, #max_time::Int64,
                0.1234 #percentage::Float64
    )
    @test tmp == [Inf, Array{Float64, 1}()]

    param_values[param_index["opt_sense"]] = MAXIMIZE
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)

    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                120, #max_time::Int64,
                0.1234 #percentage::Float64
    )
    @test tmp == [-Inf, Array{Float64, 1}()]

    ########################
    # Test maximum time
    ########################

    param_values = copy(default_param_values)
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
                5, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test ceil(time() - start_time) ≈ 6.0

    start_time = time()
    tmp = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                1, #block_size::Int64,
                10, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test ceil(time() - start_time) ≈ 11.0

    ########################
    # Test the relink
    ########################
    # **NOTE:** these tests may fail with the random number generation changes.
    # In such case, we have to figure out how to make this test better.

    data_path = joinpath(@__DIR__, "brkga_data_files")

    load_brkga_data(joinpath(data_path, "data_path_relink.jld"), brkga_data)
    results = load(joinpath(data_path, "best_solutions_pr_direct.jld"))
    brkga_data.decode! = decode!
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
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block1[1] ≈ results["block1"][1]
    @test block1[2] ≈ results["block1"][2]

    # Size 10
    chr1, chr2 = next_pair(chr2)
    block10 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block10[1] ≈ results["block10"][1]
    @test block10[2] ≈ results["block10"][2]

    # Size 100
    chr1, chr2 = next_pair(chr2)
    block100 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                100, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block100[1] ≈ results["block100"][1]
    @test block100[2] ≈ results["block100"][2]

    # Size 400
    chr1, chr2 = next_pair(chr2)
    block400 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                400, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block400[1] ≈ results["block400"][1]
    @test block400[2] ≈ results["block400"][2]

    # Size 372
    chr1, chr2 = next_pair(chr2)
    block372 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                372, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block372[1] ≈ results["block372"][1]
    @test block372[2] ≈ results["block372"][2]

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
                120, #max_time::Int64,
                0.1 #percentage::Float64
    )
    @test path10[1] ≈ results["path10"][1]
    @test path10[2] ≈ results["path10"][2]

    # Path 30%
    chr1, chr2 = next_pair(chr2)
    path30 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.3 #percentage::Float64
    )
    @test path30[1] ≈ results["path30"][1]
    @test path30[2] ≈ results["path30"][2]

    # Path 50%
    chr1, chr2 = next_pair(chr2)
    path50 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.50 #percentage::Float64
    )
    @test path50[1] ≈ results["path50"][1]
    @test path50[2] ≈ results["path50"][2]

    # Path 100%
    chr1, chr2 = next_pair(chr2)
    path100 = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test path100[1] ≈ results["path100"][1]
    @test path100[2] ≈ results["path100"][2]

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
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test xy[1] ≈ results["xy"][1]
    @test xy[2] ≈ results["xy"][2]

    # x > y
    chr1, chr2 = next_pair(chr2)
    yx = BrkgaMpIpr.direct_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> x[1] > y[2], #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test yx[1] ≈ results["yx"][1]
    @test yx[2] ≈ results["yx"][2]
end

################################################################################

@testset "permutation_based_path_relink!()" begin
    param_values = copy(default_param_values)
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
                120, #max_time::Int64,
                0.1234 #percentage::Float64
    )
    @test tmp == [-Inf, Array{Float64, 1}()]

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
                120, #max_time::Int64,
                0.1234 #percentage::Float64
    )
    @test tmp == [Inf, Array{Float64, 1}()]

    ########################
    # Test maximum time
    ########################

    param_values = copy(default_param_values)
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
                5, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test ceil(time() - start_time) ≈ 6.0

    start_time = time()
    tmp = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[1], #chromosome1
                brkga_data.current[1].chromosomes[2], #chromosome2
                (x, y) -> false, #affect_solution::Function,
                1, #block_size::Int64,
                10, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test ceil(time() - start_time) ≈ 11.0

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
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )

    # Size 10
    chr1, chr2 = next_pair(chr2)
    block10 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block10[1] ≈ results["block10"][1]
    @test block10[2] ≈ results["block10"][2]

    # Size 100
    chr1, chr2 = next_pair(chr2)
    block100 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                100, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block100[1] ≈ results["block100"][1]
    @test block100[2] ≈ results["block100"][2]

    # Size 400
    chr1, chr2 = next_pair(chr2)
    block400 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                400, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block400[1] ≈ results["block400"][1]
    @test block400[2] ≈ results["block400"][2]

    # Size 372
    chr1, chr2 = next_pair(chr2)
    block372 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                372, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test block372[1] ≈ results["block372"][1]
    @test block372[2] ≈ results["block372"][2]

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
                120, #max_time::Int64,
                0.1 #percentage::Float64
    )
    @test path10[1] ≈ results["path10"][1]
    @test path10[2] ≈ results["path10"][2]

    # Path 30%
    chr1, chr2 = next_pair(chr2)
    path30 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.3 #percentage::Float64
    )
    @test path30[1] ≈ results["path30"][1]
    @test path30[2] ≈ results["path30"][2]

    # Path 50%
    chr1, chr2 = next_pair(chr2)
    path50 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                0.5 #percentage::Float64
    )
    @test path50[1] ≈ results["path50"][1]
    @test path50[2] ≈ results["path50"][2]

    # Path 100%
    chr1, chr2 = next_pair(chr2)
    path100 = BrkgaMpIpr.permutation_based_path_relink!(
                brkga_data, #brkga_data::BrkgaData,
                brkga_data.current[1].chromosomes[chr1], #chromosome1
                brkga_data.current[1].chromosomes[chr2], #chromosome2
                (x, y) -> true, #affect_solution::Function,
                10, #block_size::Int64,
                120, #max_time::Int64,
                1.0 #percentage::Float64
    )
    @test path100[1] ≈ results["path100"][1]
    @test path100[2] ≈ results["path100"][2]

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
            120, #max_time::Int64,
            0.5 #percentage::Float64
    )
    @test xy[1] ≈ results["xy"][1]
    @test xy[2] ≈ results["xy"][2]

    # x > y
    chr1, chr2 = next_pair(chr2)
    yx = BrkgaMpIpr.permutation_based_path_relink!(
            brkga_data, #brkga_data::BrkgaData,
            brkga_data.current[1].chromosomes[chr1], #chromosome1
            brkga_data.current[1].chromosomes[chr2], #chromosome2
            (x, y) -> x[1] > y[2], #affect_solution::Function,
            10, #block_size::Int64,
            120, #max_time::Int64,
            0.5 #percentage::Float64
    )
    @test yx[1] ≈ results["yx"][1]
    @test yx[2] ≈ results["yx"][2]
end

###############################################################################

@testset "path_relink!()" begin
    param_values = copy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["num_independent_populations"]] = 1
    brkga_data = build_brkga(param_values...)

    ########################
    # Test arguments
    ########################

    @test_throws ErrorException path_relink!(brkga_data, #::BrkgaData,
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            1, #block_size::Int64,
                                            120, #max_time::Int64,
                                            1.0 #percentage::Float64
                                            )

    initialize!(brkga_data)

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            1, #block_size::Int64,
                                            120, #max_time::Int64,
                                            -0.1 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            1, #block_size::Int64,
                                            120, #max_time::Int64,
                                            1.01 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            0, #block_size::Int64,
                                            120, #max_time::Int64,
                                            1.0 #percentage::Float64
                                            )

    @test_throws ArgumentError path_relink!(brkga_data, #::BrkgaData,
                                            (x, y) -> 1.0, #compute_distance::Function,
                                            (x, y) -> true, #affect_solution::Function,
                                            10, #number_pairs::Int64,
                                            0.5, #minimum_distance::Float64,
                                            DIRECT, #::PathRelinkingType,
                                            BESTSOLUTION, # PathRelinkingSelection
                                            -10, #block_size::Int64,
                                            120, #max_time::Int64,
                                            1.0 #percentage::Float64
                                            )

    ########################
    # Test fake homogeneity
    ########################

    @test TOO_HOMOGENEOUS == path_relink!(
        brkga_data, #::BrkgaData,
        (x, y) -> 0.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        10, #number_pairs::Int64,
        1.0, #minimum_distance::Float64,
        DIRECT, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        1, #block_size::Int64,
        2, #max_time::Int64,
        1.0 #percentage::Float64
    )

    #######################
    # Test the direct path relink
    #######################

    # Let's save the data to use after this.
    original_brkga_data = deepcopy(brkga_data)

    @test BEST_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        (x, y) -> 1.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        10, #number_pairs::Int64,
        0.5, #minimum_distance::Float64,
        DIRECT, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        1, #block_size::Int64,
        10, #max_time::Int64,
        1.0, #percentage::Float64
    )

    # Do not change the best.
    @test get_best_fitness(brkga_data) > get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) != get_best_chromosome(original_brkga_data)

    ########################
    # Test the permutation path relink
    ########################

    ########################
    # No improvement found
    ########################

    param_values = copy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["num_independent_populations"]] = 1
    param_values[param_index["decode!"]] = rank_decode!
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)
    original_brkga_data = deepcopy(brkga_data)

    @test NO_IMPROVEMENT == path_relink!(
            brkga_data, #::BrkgaData,
            (x, y) -> 1.0, #compute_distance::Function,
            (x, y) -> true, #affect_solution::Function,
            10, #number_pairs::Int64,
            0.5, #minimum_distance::Float64,
            PERMUTATION, #::PathRelinkingType,
            BESTSOLUTION, # PathRelinkingSelection
            1, #block_size::Int64,
            10, #max_time::Int64,
            1.0, #percentage::Float64
        )

    # Do not change at all.
    @test get_best_fitness(brkga_data) == get_best_fitness(original_brkga_data)
    @test get_current_population(brkga_data, 1) !=
          get_current_population(original_brkga_data, 1)

    ########################
    # Improvement found, but not in the best solution.
    ########################

    param_values = copy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["pop_size"]] = 10
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 100
    param_values[param_index["instance"]] = Instance(100)
    param_values[param_index["num_independent_populations"]] = 1
    param_values[param_index["decode!"]] = rank_decode!
    brkga_data = build_brkga(param_values...)
    initialize!(brkga_data)
    original_brkga_data = deepcopy(brkga_data)

    @test ELITE_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        (x, y) -> 1.0, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        0, #number_pairs::Int64,
        0.5, #minimum_distance::Float64,
        PERMUTATION, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        1, #block_size::Int64,
        0, #max_time::Int64,
        1.0, #percentage::Float64
    )

    # Do not change the best.
    @test get_best_fitness(brkga_data) ≈ get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) ≈ get_best_chromosome(original_brkga_data)

    # But population changes.
    @test get_current_population(brkga_data, 1) !=
          get_current_population(original_brkga_data, 1)

    ########################
    # Improvement found
    ########################

    param_values = copy(default_param_values)
    param_values[param_index["seed"]] = 2700001
    param_values[param_index["pop_size"]] = 10
    param_values[param_index["opt_sense"]] = MAXIMIZE
    param_values[param_index["chr_size"]] = 10
    param_values[param_index["instance"]] = Instance(10)
    param_values[param_index["num_independent_populations"]] = 1
    param_values[param_index["decode!"]] = rank_decode!
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
    for i in 1:brkga_data.population_size
        value = brkga_data.decode!(population.chromosomes[i],
                                   brkga_data.problem_instance, false)
        population.fitness[i] = (value, i)
    end

    function dist(x, y)
        total = 0.0
        for i in 1:length(x)
            if x[i] != y[i]
                total += 1.0
            end
        end
        return total
    end

    original_brkga_data = deepcopy(brkga_data)

    @test BEST_IMPROVEMENT == path_relink!(
        brkga_data, #::BrkgaData,
        dist, #compute_distance::Function,
        (x, y) -> true, #affect_solution::Function,
        0, #number_pairs::Int64,
        1.0, #minimum_distance::Float64,
        PERMUTATION, #::PathRelinkingType,
        RANDOMELITE, # PathRelinkingSelection
        1, #block_size::Int64,
        0, #max_time::Int64,
        1.0, #percentage::Float64
    )

    @test get_best_fitness(brkga_data) > get_best_fitness(original_brkga_data)
    @test get_best_chromosome(brkga_data) != get_best_chromosome(original_brkga_data)
end
