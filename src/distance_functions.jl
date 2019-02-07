################################################################################
# distance_functions.jl: Definitions of functions that compute the distance
# between two vectors of double numbers.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
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

################################################################################
# Hamming distance functions
################################################################################

"""
    hamming_distance(vector1::Array{Float64, 1}, vector2::Array{Float64, 1},
                     threshold::Float64 = 0.5)::Float64

Compute the [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance)
between two vectors. It takes a `threshold` parameter to "binarize" the
vectors. For instance, if `threshold = 0.7`, all values larger than or equal
to 0.7 will be considerd `1.0`, otherwise `0.0`.

**NOTE:** this function may be more appropriated to threshold/direct
chromosome representations.

# Arguments
- `vector1::Array{Float64, 1}`: the first vector.
- `vector2::Array{Float64, 1}`: the second vector.
- `threshold::Float64 = 0.5`: the threshold for binarization.

# Throws
- `ArgumentError`: if `vector1` and `vector2` have different sizes.
"""
function hamming_distance(vector1::Array{Float64, 1},
                          vector2::Array{Float64, 1};
                          threshold::Float64 = 0.5)::Float64
    if length(vector1) != length(vector2)
        throw(ArgumentError("Input vectors have different sizes."))
    end

    dist::Int64 = 0
    @inbounds for i in 1:length(vector1)
        if (vector1[i] < threshold) != (vector2[i] < threshold)
            dist += 1
        end
    end
    return Float64(dist)
end

################################################################################

"""
    affect_solution_hamming_distance(block1::SubArray{Float64, 1},
                                     block2::SubArray{Float64, 1},
                                     threshold::Float64 = 0.5)::Bool

Return `true` the the changing of the blocks of keys `block1` by the
blocks of keys `block2` affects the solution, based on the
[Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance).

**NOTE 1:** this function may be more appropriated to threshold/direct
chromosome representations.

**NOTE 2:** `block1` and `block2` must have the same size. No bounds checking
is done due to performance reasons.

**NOTE 3:** this function is annotated with `@inline` due to performance
reasons too.

# Arguments
- `block1::SubArray{Float64, 1}`: the first vector.
- `block2::SubArray{Float64, 1}`: the second vector.
- `threshold::Float64 = 0.5`: the threshold for binarization.
"""
@inline function affect_solution_hamming_distance(
            block1::SubArray{Float64, 1},
            block2::SubArray{Float64, 1};
            threshold::Float64 = 0.5)::Bool
    for i in 1:length(block1)
        if (block1[i] < threshold) != (block2[i] < threshold)
            return true
        end
    end
    return false
end

################################################################################
# Kendall Tau distance functions
################################################################################

"""
    kendall_tau_distance(vector1::Array{Float64, 1},
                         vector2::Array{Float64, 1};
                         stop_immediately::Bool = false)::Float64

    kendall_tau_distance(vector1::SubArray{Float64, 1},
                         vector2::SubArray{Float64, 1};
                         stop_immediately::Bool = false)::Float64

Compute the [Kendall Tau distance](https://en.wikipedia.org/wiki/Kendall_tau_distance)
between two vectors.

**NOTE:** this function may be more appropriated to permutation chromosome
representations.

# Arguments
- `vector1::Array{Float64, 1}`: the first vector.
- `vector2::Array{Float64, 1}`: the second vector.
- `stop_immediately::Bool = false`: if `true`, stop the computation
  immediately after find a difference.

# Throws
- `ArgumentError`: if `vector1` and `vector2` have different sizes.
"""
function kendall_tau_distance(vector1::Array{Float64, 1},
                              vector2::Array{Float64, 1};
                              stop_immediately::Bool = false)::Float64
    return kendall_tau_distance(view(vector1, 1:length(vector1)),
                                view(vector2, 1:length(vector2)),
                                stop_immediately = stop_immediately)
end

function kendall_tau_distance(vector1::SubArray{Float64, 1},
                              vector2::SubArray{Float64, 1};
                              stop_immediately::Bool = false)::Float64
    if length(vector1) != length(vector2)
        throw(ArgumentError("Input vectors have different sizes."))
    end

    pairs_v1 = Array{Tuple{Float64, Int64}}(undef, length(vector1))
    @inbounds for (rank, key) in enumerate(vector1)
        pairs_v1[rank] = (key, rank)
    end
    sort!(pairs_v1)

    pairs_v2 = Array{Tuple{Float64, Int64}}(undef, length(vector2))
    @inbounds for (rank, key) in enumerate(vector2)
        pairs_v2[rank] = (key, rank)
    end
    sort!(pairs_v2)

    disagreements::Int64 = 0
    for i in 1:(length(vector1) - 1), j in (i + 1):length(vector1)
        if ((pairs_v1[i][2] < pairs_v1[j][2]
            && pairs_v2[i][2] > pairs_v2[j][2])) ||
           (pairs_v1[i][2] > pairs_v1[j][2] && pairs_v2[i][2] < pairs_v2[j][2])

            disagreements += 1
            if stop_immediately
                break
            end
        end
    end

    return Float64(disagreements)
end

################################################################################

"""
    affect_solution_kendall_tau(block1::SubArray{Float64, 1},
                                block2::SubArray{Float64, 1})::Bool

Return `true` the the changing of the blocks of keys `block1` by the
blocks of keys `block2` affects the solution, based on the
[Kendall Tau distance.](https://en.wikipedia.org/wiki/Kendall_tau_distance)

**NOTE:** `block1` and `block2` must have the same size. No bounds checking
is done due to performance reasons.

# Arguments
- `block1::SubArray{Float64, 1}`: the first vector.
- `block2::SubArray{Float64, 1}`: the second vector.
"""
function affect_solution_kendall_tau(
            block1::SubArray{Float64, 1},
            block2::SubArray{Float64, 1})::Bool

    if length(block1) == 1
        return abs(block1[1] - block2[1]) > 1e-6
    end

    return kendall_tau_distance(block1, block2, stop_immediately = true) > 1e-6
end
