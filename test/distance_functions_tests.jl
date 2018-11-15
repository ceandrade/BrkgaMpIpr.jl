################################################################################
# distance_functions_tests.jl.jl: unit tests for distance functions.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Nov 15, 2018 by ceandrade
# Last update: Nov 15, 2018 by ceandrade
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

@testset "hamming_distance()" begin
    # Create and warm up a random number generator.
    local_rng = MersenneTwister(2700001)
    rand(local_rng, 1000)

    # Exceptions
    x = rand(local_rng, 10)
    y = rand(local_rng, 5)

    @test_throws ArgumentError hamming_distance(x, y)

    # Manual test 1
    x = [0.3, 0.4, 0.7]
    y = [0.4, 0.6, 0.7]

    @test hamming_distance(x, y) ≈ 1.0
    @test hamming_distance(x, y, threshold = 0.1) ≈ 0.0
    @test hamming_distance(x, y, threshold = 0.9) ≈ 0.0

    # Manual test 2
    x = [0.1, 0.2, 0.3]
    y = [0.6, 0.7, 0.8]

    @test hamming_distance(x, y) ≈ 3.0
    @test hamming_distance(x, y, threshold = 0.1) ≈ 0.0
    @test hamming_distance(x, y, threshold = 0.9) ≈ 0.0
    @test hamming_distance(x, y, threshold = 0.65) ≈ 2.0

    # Random numbers
    x = rand(local_rng, 1000)
    y = rand(local_rng, 1000)

    @test hamming_distance(x, y, threshold = 0.5) ≈ 547.0
    @test hamming_distance(x, y, threshold = 0.25) ≈ 384.0
    @test hamming_distance(x, y, threshold = 0.75) ≈ 364.0
    @test hamming_distance(x, y, threshold = 0.10) ≈ 177.0
    @test hamming_distance(x, y, threshold = 0.90) ≈ 171.0
end

################################################################################

@testset "affect_solution_hamming_distance()" begin
    x = [0.1, 0.2, 0.4, 0.6, 0.7]
    y = [0.1, 0.3, 0.5, 0.6, 0.7]

    rg = 1:1
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == false

    rg = 3:3
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == true

    @test affect_solution_hamming_distance(view(x, rg), view(y, rg),
                                           threshold = 0.6) == false

    rg = 2:2
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == false

    rg = 1:2
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == false

    rg = 2:3
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == true

    rg = 3:4
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == true

    rg = 4:5
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == false

    rg = 1:length(x)
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == true

    rg = 1:length(x)
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg),
                                           threshold = 0.1) == false

    @test affect_solution_hamming_distance(view(x, rg), view(y, rg),
                                           threshold = 0.9) == false
end

################################################################################

@testset "kendall_tau_distance()" begin
    # Create and warm up a random number generator.
    local_rng = MersenneTwister(2700001)
    rand(local_rng, 1000)

    # Exceptions
    x = rand(local_rng, 10)
    y = rand(local_rng, 5)

    @test_throws ArgumentError kendall_tau_distance(x, y)

    # Manual test 1
    x = [0.3, 0.4, 0.7]
    y = [0.1, 0.2, 0.3]
    @test kendall_tau_distance(x, y) ≈ 0.0

    # Manual test 2
    x = [0.3, 0.4, 0.7]
    y = [0.7, 0.4, 0.6]
    @test kendall_tau_distance(x, y) ≈ 2.0

    # Random test 1
    x = rand(local_rng, 1000)
    y = rand(local_rng, 1000)
    @test kendall_tau_distance(x, y) ≈ 250058.0

    # Random test 2
    x = sort(rand(local_rng, 1000))
    y = sort(rand(local_rng, 1000))
    @test kendall_tau_distance(x, y) ≈ 0.0

    # Random test 3 (all pairs have difference)
    x = sort(rand(local_rng, 1000))
    y = sort(rand(local_rng, 1000), rev = true)
    @test kendall_tau_distance(x, y) ≈ (1000 * 999) / 2.0
end

################################################################################

@testset "affect_solution_kendall_tau()" begin
    x = [0.1, 0.2, 0.4, 0.6, 0.7]
    y = [0.1, 0.3, 0.5, 0.6, 0.7]

    rg = 1:1
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == false

    rg = 3:3
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == true

    rg = 5:5
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == false

    rg = 1:2
    @test affect_solution_hamming_distance(view(x, rg), view(y, rg)) == false

    rg = 2:4
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == false

    x = [0.1, 0.2, 0.4, 0.6, 0.7]
    y = [0.1, 0.3, 0.5, 0.3, 0.1]

    rg = 1:3
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == false

    rg = 1:4
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == true

    # Create and warm up a random number generator.
    local_rng = MersenneTwister(2700001)
    rand(local_rng, 1000)
    x = rand(local_rng, 1000)
    y = rand(local_rng, 1000)

    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == true

    x = sort(rand(local_rng, 1000))
    y = sort(rand(local_rng, 1000))
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == false

    x = sort(rand(local_rng, 1000))
    y = sort(rand(local_rng, 1000), rev = true)
    @test affect_solution_kendall_tau(view(x, rg), view(y, rg)) == true
end
