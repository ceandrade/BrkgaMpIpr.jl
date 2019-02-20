################################################################################
# path_relink.jl: main evolutionary routines.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 06, 2018 by ceandrade
# Last update: Feb 19, 2019 by ceandrade
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

import Base

################################################################################
# Auxiliary structures and functions
################################################################################

"""
    mutable struct Triple

Hold the data structures used to build a candidate chromosome for parallel
decoding on direct path relink.

!!! warning
    THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.
"""
mutable struct Triple
    chr::Array{Float64, 1}
    fitness::Float64
    block_index::Int64
    Triple() = new(Array{Float64, 1}(), 0.0, 0)
    Triple(chr_::Array{Float64, 1}) = new(chr_, 0.0, 0)
end

################################################################################

"""
    mutable struct DecodeStruct

Hold the data structures used to build a candidate chromosome for parallel
decoding on permutation-based path relink.

!!! warning
    THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.
"""
mutable struct DecodeStruct
    chr::Array{Float64, 1}
    fitness::Float64
    key_index::Int64
    pos1::Int64
    pos2::Int64
    DecodeStruct() = new(Array{Float64, 1}(), 0.0, 0, 0, 0)
    DecodeStruct(chr_::Array{Float64, 1}) = new(chr_, 0.0, 0, 0, 0)
end

################################################################################

"""
    find_block_range(block_number::Int64, block_size::Int64,
                     max_end::Int64)::UnitRange{Int64}

Return a positive range for the given `block_number` with length `block_size`,
limited to the `max_end`.

!!! note
    This function only accept positive numbers, and all sanity check is
    disregarded due to performance reasons.

!!! warning
    THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.
"""
@inline function find_block_range(block_number::Int64, block_size::Int64,
                                  max_end::Int64)::UnitRange{Int64}
     block_base = block_size * (block_number - 1) + 1
     block_end = min((block_base + block_size - 1), max_end)
    return block_base:block_end
end

################################################################################

"""
    swap!(x::Array{Any, 1}, pos1::Int64, pos2::Int64)

Swap the value in position `pos1` with the value in position `pos2` in
vector `x`.

!!! note
    This function only accept positive numbers, and all sanity and bounds check
    is disregarded due to performance reasons.

!!! warning
    THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.
"""
@inline function swap!(x::Array{T, 1}, pos1::Int64, pos2::Int64) where T
    @inbounds x[pos1], x[pos2] = x[pos2], x[pos1]
    nothing
end

################################################################################

"""
    Base.:|(x::PathRelinkingResult,
            y::PathRelinkingResult)::PathRelinkingResult

Perform bitwise `OR` between two [`PathRelinkingResult`](@ref) returning
the highest rank `PathRelinkingResult`.

# Examples

```julia-repl
julia> TOO_HOMOGENEOUS | NO_IMPROVEMENT
NO_IMPROVEMENT::PathRelinkingResult = 1

julia> NO_IMPROVEMENT | ELITE_IMPROVEMENT
ELITE_IMPROVEMENT::PathRelinkingResult = 3

julia> ELITE_IMPROVEMENT | BEST_IMPROVEMENT
BEST_IMPROVEMENT::PathRelinkingResult = 7
```
"""
@inline function Base.:|(x::PathRelinkingResult,
                         y::PathRelinkingResult)::PathRelinkingResult
    return PathRelinkingResult(Int64(x) | Int64(y))
end

################################################################################
# Direct path relinking method
################################################################################

"""
    function direct_path_relink!(brkga_data::BrkgaData,
                                 chromosome1::Array{Float64, 1},
                                 chromosome2::Array{Float64, 1},
                                 affect_solution::Function,
                                 block_size::Int64,
                                 max_time::Float64,
                                 percentage::Float64
        )::Tuple{Float64, Array{Float64, 1}}

Perform the direct path relinking, changing each allele or block
of alleles of base chromosome for the correspondent one in the guide
chromosome.

The API will call `decode!()` function always with `writeback = false`. The
reason is that if the decoder rewrites the chromosome, the path between
solutions is lost and inadvertent results may come up. Note that at the end
of the path relinking, the method calls the decoder with `writeback = true`
in the best chromosome found to guarantee that this chromosome is re-written
to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

!!! warning
    As it is in [`evolve!()`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

!!! warning
    THIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS
    CALLED FROM THE [`path_relink!()`](@ref) FUNCTION. Due to this reason,
    this method **DOES NOT** perform health checks on the arguments.

# Arguments
- [`brkga_data::BrkgaData`](@ref BrkgaData): the BRKGA data.

- `chromosome1::Array{Float64, 1}` and `chromosome2::Array{Float64, 1}`: the
  chromosomes to be used to build the path.

- `affect_solution::Function`: function that takes two partial chromosomes /
  block of genes `block1` and `block2` and checks whether changing the keys from
  `block1` to `block2` affects the solution. For instance, suppose that the
  alleles/keys are used as threshold such that values > 0.5 activate a feature.
  Suppose we have `block1 = [0.3, 0.4, 0.1]` and `block2 = [0.4, 0.1, 0.2]`.
  Since all values are below 0.5, changing the keys from `block1` to `block2`
  does not chage the solution, and therefore, we can drop such change (and
  subsequentely decoding). The blocks can hold only one key/allele, sequential
  key blocks, of even the whole chromosome. `affect_solution` takes two
  views/subarrays. The function **MUST HAVE** the following signature

  ```julia
  affect_solution(block1::SubArray{Float64, 1},
                  block2::SubArray{Float64, 1})::Bool
  ```

  !!! note
        This function depends on the problem structure and how the
        keys/alleles are used.

- `block_size::Int64`: (posite) number of alleles to be exchanged at once in
  each iteration. If `block_size == 1`, the traditional path relinking is
  performed.

- `max_time::Float64`: abort path-relinking when reach `max_time`.
  If `max_time <= 0`, no limit is imposed. Given in seconds.

- `percentage::Float64`: define the size, in percentage, of the path to
  build. Range [0, 1].

# Returns

- `Array{Any, 1}`: the best pair [fitness, chromosome]
  found during the relinking. If the relink is not possible due to homogeneity,
  `-Inf` returns in case of maximization, and `Inf` in case of minimization.

"""
function direct_path_relink!(brkga_data::BrkgaData,
                             chromosome1::Array{Float64, 1},
                             chromosome2::Array{Float64, 1},
                             affect_solution::Function,
                             block_size::Int64,
                             max_time::Float64,
                             percentage::Float64
    )::Array{Any, 1}

    pr_start_time = time()
    bd = brkga_data

    best_chr_found = Array{Float64, 1}(undef, bd.chromosome_size)
    best_fitness_found::Float64 = (bd.opt_sense == MAXIMIZE) ? -Inf : Inf

    num_blocks = cld(bd.chromosome_size, block_size)
    path_size = Int64(floor(percentage * num_blocks))
    remaining_blocks = collect(1:num_blocks)

    # Allocate memory for the candidates.
    candidates_base = Array{Triple, 1}(undef, num_blocks)
    candidates_guide = Array{Triple, 1}(undef, num_blocks)

    # References for the base and guide chromosomes.
    base = chromosome1
    guide = chromosome2

    # **NOTE:** two loops are faster than one according to `@benchmark`,
    # probably because of cache racing.
    @inbounds Threads.@threads for i in 1:num_blocks
        candidates_base[i] = Triple(copy(base))
    end
    @inbounds Threads.@threads for i in 1:num_blocks
        candidates_guide[i] = Triple(copy(guide))
    end

    # Holds the original keys.
    old_keys = Array{Float64, 1}(undef, bd.chromosome_size)

    sense = (bd.opt_sense == MAXIMIZE)
    iterations = 1
    while !isempty(remaining_blocks)
        state = iterate(remaining_blocks)

        i = 1
        while state !== nothing && !isempty(remaining_blocks)
            it_block_idx = state[2] - 1
            block_range = find_block_range(remaining_blocks[it_block_idx],
                                           block_size, length(guide))
            block1 = view(candidates_base[i].chr, block_range)
            block2 = view(guide, block_range)

            if !affect_solution(block1, block2)
                deleteat!(remaining_blocks, it_block_idx)
                if it_block_idx > length(remaining_blocks)
                    state = nothing
                end
                continue
            end

            # Save the former keys before...
            old_keys[block_range] = candidates_base[i].chr[block_range]

            # ... copy the keys from the guide solution.
            candidates_base[i].chr[block_range] = guide[block_range]

            candidates_base[i].block_index = it_block_idx
            state = iterate(remaining_blocks, state[2])
            i += 1
        end

        if i == 1 && isempty(remaining_blocks)
            break
        end

        Threads.@threads for i in 1:length(remaining_blocks)
            candidates_base[i].fitness = (sense) ? -Inf : Inf
            if time() - pr_start_time > max_time
                continue
            end

            candidates_base[i].fitness =
                bd.decode!(candidates_base[i].chr, bd.problem_instance, false)
        end

        # Locate the best candidate.
        best_value::Float64 = sense ? -Inf : Inf
        best_index = 1
        best_block_index = 1

        @inbounds for i in 1:length(remaining_blocks)
            if ((best_value < candidates_base[i].fitness && sense) ||
                (best_value > candidates_base[i].fitness && !sense))
                best_value = candidates_base[i].fitness
                best_index = i
                best_block_index = candidates_base[i].block_index
            end
        end

        # Hold it, if it is the best found until now.
        if ((sense && best_fitness_found < candidates_base[best_index].fitness)
            ||
            (!sense && best_fitness_found > candidates_base[best_index].fitness))
            best_fitness_found = candidates_base[best_index].fitness
            best_chr_found[:] = candidates_base[best_index].chr
        end

        # Restore original keys and copy the block of keys for all future
        # candidates. The last candidate will not be used.
        # it_block_idx = start(remaining_blocks)
        for i in 1:(length(remaining_blocks) - 1)
            it_block_idx = iterate(remaining_blocks)[2] - 1
            block_range = find_block_range(remaining_blocks[it_block_idx],
                                           block_size, length(guide))

            candidates_base[i].chr[block_range] = old_keys[block_range]

            block_range = find_block_range(best_block_index, block_size,
                                           length(guide))
            candidates_base[i].chr[block_range] =
                candidates_base[best_index].chr[block_range]
        end

        # Swap roles.
        base, guide = guide, base
        candidates_base, candidates_guide = candidates_guide, candidates_base
        deleteat!(remaining_blocks, best_block_index)

        iterations += 1
        if iterations == path_size || time() - pr_start_time > max_time
            break
        end
    end

    if best_fitness_found == Inf || best_fitness_found == -Inf
        best_chr_found = Array{Float64, 1}()
    end
    return [best_fitness_found, best_chr_found]
end

################################################################################
# Permutation-based path relinking method
################################################################################

"""
    function permutation_based_path_relink!(brkga_data::BrkgaData,
                                            chromosome1::Array{Float64, 1},
                                            chromosome2::Array{Float64, 1},
                                            affect_solution::Function,
                                            block_size::Int64,
                                            max_time::Float64,
                                            percentage::Float64
        )::Tuple{Float64, Array{Float64, 1}}

Perform the permutation-based path relinking. In this method, the permutation
induced by the keys in the guide solution is used to change the order of the
keys in the permutation induced by the base solution.

The API will call `decode!()` function always with `writeback = false`. The
reason is that if the decoder rewrites the chromosome, the path between
solutions is lost and inadvertent results may come up. Note that at the end
of the path relinking, the method calls the decoder with `writeback = true`
in the best chromosome found to guarantee that this chromosome is re-written
to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

!!! warning
    As it is in [`evolve!()`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

!!! warning
    THIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS
    CALLED FROM THE [`path_relink!()`](@ref) FUNCTION. Due to this reason,
    this method **DOES NOT** perform health checks on the arguments.

# Arguments
- [`brkga_data::BrkgaData`](@ref BrkgaData): the BRKGA data.

- `chromosome1::Array{Float64, 1}` and `chromosome2::Array{Float64, 1}`: the
  chromosomes to be used to build the path.

- `affect_solution::Function`: not used in this function but kept to API
  compatibility.

- `block_size::Int64`: not used in this function but kept to API compatibility.

- `max_time::Float64`: abort path-relinking when reach `max_time`.
  If `max_time <= 0`, no limit is imposed. Given in seconds.

- `percentage::Float64`: define the size, in percentage, of the path to
       build. Range [0, 1].

# Returns

- `Array{Any, 1}`: the best pair [fitness, chromosome]
  found during the relinking. If the relink is not possible due to homogeneity,
  `-Inf` returns in case of maximization, and `Inf` in case of minimization.
"""
function permutation_based_path_relink!(brkga_data::BrkgaData,
                                        chromosome1::Array{Float64, 1},
                                        chromosome2::Array{Float64, 1},
                                        affect_solution::Function,
                                        block_size::Int64,
                                        max_time::Float64,
                                        percentage::Float64
    )::Array{Any, 1}

    bd = brkga_data
    pr_start_time = time()

    best_chr_found = Array{Float64, 1}(undef, bd.chromosome_size)
    best_fitness_found::Float64 = (bd.opt_sense == MAXIMIZE) ? -Inf : Inf

    path_size = Int64(floor(percentage * bd.chromosome_size))
    remaining_indices = collect(1:bd.chromosome_size)

    # Allocate memory for the candidates.
    candidates_base = Array{DecodeStruct, 1}(undef, bd.chromosome_size)
    candidates_guide = Array{DecodeStruct, 1}(undef, bd.chromosome_size)

    # References for the base and guide chromosomes.
    base = chromosome1
    guide = chromosome2

    # **NOTE:** two loops are faster than one according to `@benchmark`,
    # probably because of cache racing.
    @inbounds Threads.@threads for i in 1:bd.chromosome_size
        candidates_base[i] = DecodeStruct(copy(base))
    end
    @inbounds Threads.@threads for i in 1:bd.chromosome_size
        candidates_guide[i] = DecodeStruct(copy(guide))
    end

    # Create and order the indices.
    base_indices = collect(1:bd.chromosome_size)
    guide_indices = collect(1:bd.chromosome_size)
    sorted = Array{Tuple{Float64, Int64}, 1}(undef, bd.chromosome_size)

    for j in 1:2
        for i in 1:bd.chromosome_size
            sorted[i] = (base[i], i)
        end
        sort!(sorted)

        for i in 1:bd.chromosome_size
            base_indices[i] = sorted[i][2]
        end

        base, guide = guide, base
        base_indices, guide_indices = guide_indices, base_indices
    end
    base, guide = guide, base
    base_indices, guide_indices = guide_indices, base_indices

    sense = (bd.opt_sense == MAXIMIZE)
    iterations = 1
    while !isempty(remaining_indices)
        state = iterate(remaining_indices)

        i = 1
        while state !== nothing && !isempty(remaining_indices)
            it_idx = state[2] - 1
            position_in_base = base_indices[remaining_indices[it_idx]]
            position_in_guide = guide_indices[remaining_indices[it_idx]]

            if position_in_base == position_in_guide
                deleteat!(remaining_indices, it_idx)
                if it_idx > length(remaining_indices)
                    state = nothing
                end
                continue
            end

            candidates_base[i].key_index = it_idx
            candidates_base[i].pos1 = position_in_base
            candidates_base[i].pos2 = position_in_guide
            candidates_base[i].fitness = (sense) ? -Inf : Inf

            state = iterate(remaining_indices, state[2])
            i += 1
        end

        if i == 1 && isempty(remaining_indices)
            break
        end

        # Decode the candidates.
        Threads.@threads for i in 1:length(remaining_indices)
            if time() - pr_start_time > max_time
                continue
            end

            swap!(candidates_base[i].chr,
                  candidates_base[i].pos1, candidates_base[i].pos2)

            candidates_base[i].fitness =
                bd.decode!(candidates_base[i].chr, bd.problem_instance, false)

            swap!(candidates_base[i].chr,
                  candidates_base[i].pos1, candidates_base[i].pos2)
        end

        # Locate the best candidate.
        best_value::Float64 = (sense) ? -Inf : Inf
        best_index = 1
        best_key_index = 1

        @inbounds for i in 1:length(remaining_indices)
            if ((best_value < candidates_base[i].fitness && sense) ||
                (best_value > candidates_base[i].fitness && !sense))
                best_value = candidates_base[i].fitness
                best_index = i
                best_key_index = candidates_base[i].key_index
            end
        end

        position_in_base = base_indices[best_key_index];
        position_in_guide = guide_indices[best_key_index];

        # Commit the best exchange in all candidates.
        # The last candidate will not be used.
        @inbounds for i in 1:length(remaining_indices) - 1
            swap!(candidates_base[i].chr, position_in_base, position_in_guide)
        end

        swap!(base_indices, position_in_base, position_in_guide)

        # Hold it, if it is the best found until now.
        if ((sense && best_fitness_found < candidates_base[best_index].fitness)
            ||
            (!sense && best_fitness_found > candidates_base[best_index].fitness))
            best_fitness_found = candidates_base[best_index].fitness
            best_chr_found .= candidates_base[best_index].chr
        end

        # Swap roles.
        base, guide = guide, base
        candidates_base, candidates_guide = candidates_guide, candidates_base
        deleteat!(remaining_indices, best_key_index)

        iterations += 1
        if iterations == path_size || time() - pr_start_time > max_time
            break
        end
    end

    if best_fitness_found == Inf || best_fitness_found == -Inf
        best_chr_found = Array{Float64, 1}()
    end
    return [best_fitness_found, best_chr_found]
end

################################################################################
# Main path relinking method
################################################################################

"""
    function path_relink!(brkga_data::BrkgaData,
                          pr_type::PathRelinkingType,
                          pr_selection::PathRelinkingSelection,
                          compute_distance::Function,
                          affect_solution::Function,
                          minimum_distance::Float64,
                          block_size::Int64,
                          max_time::Int64,
                          percentage::Float64
    )::PathRelinkingResult

Perform path relinking between elite solutions that are, at least, a given
minimum distance between themselves. In this method, the local/loaded
parameters are ignored in favor to the supplied ones.

In the presence of multiple populations, the path relinking is performed
between elite chromosomes from different populations, in a circular fashion.
For example, suppose we have 3 populations. The framework performs 3 path
relinkings: the first between individuals from populations 1 and 2, the
second between populations 2 and 3, and the third between populations 3 and 1.
In the case of just one population, both base and guiding individuals are
sampled from the elite set of that population.

Note that the algorithm tries to find a pair of base and guiding solutions
with a minimum distance given by the distance function. If this is not
possible, a new pair of solutions are sampled (without replacement) and
tested against the distance. In case it is not possible to find such pairs
for the given populations, the algorithm skips to the next pair of
populations (in a circular fashion, as described above). Yet, if such pairs
are not found in any case, the algorithm declares failure. This indicates
that the populations are very homogeneous.

If the found solution is the best solution found so far, IPR replaces the
worst solution by it. Otherwise, IPR computes the distance between the found
solution and all other solutions in the elite set, and replaces the worst
solution by it if and only if the found solution is, at least,
`minimum_distance` from all them.

The API will call `decode!()` function always with `writeback = false`. The
reason is that if the decoder rewrites the chromosome, the path between
solutions is lost and inadvertent results may come up. Note that at the end
of the path relinking, the method calls the decoder with `writeback = true`
in the best chromosome found to guarantee that this chromosome is re-written
to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

!!! warning
    As it is in [`evolve!()`](@ref), the decoding is done in parallel using
    threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
    If such property cannot be held, we suggest using single thread by
    setting the environmental variable `JULIA_NUM_THREADS = 1` [(see Julia
    Parallel Computing)]
    (https://docs.julialang.org/en/v1.1/manual/parallel-computing/).

# Arguments
- [`brkga_data::BrkgaData`](@ref BrkgaData): the BRKGA data.

- [`pr_type::PathRelinkingType`](@ref PathRelinkingType): type of path relinking
  to be performed. Either `DIRECT` or `PERMUTATION`-based.

- [`pr_selection::PathRelinkingSelection`](@ref PathRelinkingSelection):
  selection of which individuals use to path relinking. Either `BESTSOLUTION`
  or `RANDOMELITE`.

- `compute_distance::Function`: the function used to compute the distance
  between two chromosomes. The function **MUST HAVE** the following signature

  ```julia
  compute_distance(vector1::Array{Float64, 1},
                   vector2::Array{Float64, 1})::Float64
  ```

- `affect_solution::Function`: function that takes two partial chromosomes /
  block of genes `block1` and `block2` and checks whether changing the keys from
  `block1` to `block2` affects the solution. For instance, suppose that the
  alleles/keys are used as threshold such that values > 0.5 activate a feature.
  Suppose we have `block1 = [0.3, 0.4, 0.1]` and `block2 = [0.4, 0.1, 0.2]`.
  Since all values are below 0.5, changing the keys from `block1` to `block2`
  do not change the solution, and therefore, we can drop such change (and
  subsequentely decoding). The blocks can hold only one key/allele, sequential
  key blocks, or even the whole chromosome. `affect_solution` takes two
  views/subarrays. The function **MUST HAVE** the following signature

  ```julia
  affect_solution(block1::SubArray{Float64, 1},
                  block2::SubArray{Float64, 1})::Bool
  ```

  !!! note
        This function depends on the problem structure and how the
        keys/alleles are used

 - `number_pairs::Int64`: number of chromosome pairs to be tested.
   If `number_pairs < 1`, all pairs are tested.

- `minimum_distance::Float64`: minimum distance between two chromosomes computed
  by `compute_distance`.

- `block_size::Int64 = 1`: number of alleles to be exchanged at once in each
  iteration. If one, the traditional path relinking is performed.
  It must be ≥ 1.

- `max_time::Float64 = 0`: abort path-relinking when reach `max_time`.
  If `max_time ≤ 0`, no limit is imposed. Given in seconds.

- `percentage::Float64 = 1.0`: define the size, in percentage, of the path to
  build. Range [0, 1].

# Returns

- Returns [`PathRelinkingResult`](@ref) depending of the relink status.

# Throws
- `ErrorException`: if [`initialize!()`](@ref) was not called before.
- `ArgumentError`: when `percentage < 1e-6 || percentage > 1.0` and
  `block_size < 1`.
"""
function path_relink!(brkga_data::BrkgaData,
                      pr_type::PathRelinkingType,
                      pr_selection::PathRelinkingSelection,
                      compute_distance::Function,
                      affect_solution::Function,
                      number_pairs::Int64,
                      minimum_distance::Float64,
                      block_size::Int64 = 1,
                      max_time::Float64 = 0.0,
                      percentage::Float64 = 1.0
    )::PathRelinkingResult

    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. " *
              "Call initialize!() before path_relink!()")
    end

    if percentage < 1e-6 || percentage > 1.0
        throw(ArgumentError("Percentage/size of path relinking invalid: " *
                            "$(percentage)"))
    end

    if block_size < 1
        throw(ArgumentError("Invalid block size: $(block_size)"))
    end

    if max_time <= 0
        max_time = typemax(Float64)
    end

    bd = brkga_data
    initial_solution = Array{Float64, 1}(undef, bd.chromosome_size)
    guiding_solution = Array{Float64, 1}(undef, bd.chromosome_size)

    # Perform path relinking between elite chromosomes from different
    # populations. This is done in a circular fashion.
    path_relinking_possible::Bool = false
    pr_start_time = time()
    pop_count = 1
    final_status = TOO_HOMOGENEOUS
    while pop_count <= bd.params.num_independent_populations
        elapsed_seconds = time() - pr_start_time
        if elapsed_seconds > max_time
            break
        end

        pop_base = pop_count
        pop_count += 1
        pop_guide = pop_count
        found_pair = false

        # If we have just one population, we take the both solution from it.
        if bd.params.num_independent_populations == 1
            pop_base = pop_guide = 1
            pop_count = Inf

        # If we have two populations, perform just one path relinking.
        elseif bd.params.num_independent_populations == 2
            pop_count = Inf
        end

        # Do the circular thing.
        if pop_guide == bd.params.num_independent_populations + 1
            pop_guide = 1
        end

        index_pairs = Array{Pair{Int64, Int64}}(undef,
                                                bd.elite_size * bd.elite_size)

        @inbounds for i in 1:bd.elite_size, j in 1:bd.elite_size
            index_pairs[(i - 1) * bd.elite_size + j] = Pair(i, j)
        end

        tested_pairs_count = 0
        if number_pairs == 0
            number_pairs = length(index_pairs)
        end

        while !isempty(index_pairs) && tested_pairs_count < number_pairs &&
              elapsed_seconds < max_time

            index = (pr_selection == BESTSOLUTION) ? 1 :
                    rand(bd.rng, 1:length(index_pairs))
            (pos1, pos2) = index_pairs[index]

            tmp = bd.current[pop_base]
            chr1 = tmp.chromosomes[tmp.fitness[pos1][2]]

            tmp = bd.current[pop_guide]
            chr2 = tmp.chromosomes[tmp.fitness[pos2][2]]

            if compute_distance(chr1, chr2) >= minimum_distance
                initial_solution .= chr1
                guiding_solution .= chr2
                found_pair = true
                break
            end

            tested_pairs_count += 1
            elapsed_seconds = time() - pr_start_time
        end

        # The elite sets are too homogeneous, we cannot do
        #  a good path relinking. Let's try other populations.
        path_relinking_possible |= found_pair
        if !found_pair
            continue
        end

        # Perform the path relinking.
        func = (pr_type == DIRECT) ? direct_path_relink! :
                                     permutation_based_path_relink!

        # best_found[1] -> fitness
        # best_found[2] -> chromosome
        best_found = func(bd, initial_solution, guiding_solution,
                          affect_solution, block_size,
                          max_time - elapsed_seconds, percentage)

        final_status |= NO_IMPROVEMENT

        if (best_found[1] == Inf || best_found[1] == -Inf) &&
            length(best_found[2]) == 0
            continue
        end

        # Re-decode and apply local search if the decoder are able to do it.
        best_found[1] = bd.decode!(best_found[2], bd.problem_instance, true)

        # Now, check if the best solution found is really good.
        # If it is the best, overwrite the worse solution in the population.
        current = bd.current[pop_base]

        sense = bd.opt_sense == MAXIMIZE
        include_in_population =
           ( sense && best_found[1] > current.fitness[1][1]) ||
           (!sense && best_found[1] < current.fitness[1][1])

        best_overall = get_best_fitness(bd)
        if ( sense && best_found[1] > best_overall) ||
           (!sense && best_found[1] < best_overall)
            final_status |= BEST_IMPROVEMENT
        end

        if (!include_in_population &&
            (( sense && best_found[1] > current.fitness[bd.elite_size][1]) ||
             (!sense && best_found[1] < current.fitness[bd.elite_size][1])))

            include_in_population = true
            @inbounds for i in 1:bd.elite_size
                if compute_distance(best_found[2],
                                    current.chromosomes[current.fitness[i][2]]
                   ) < minimum_distance - 1e-6
                    include_in_population = false
                    final_status |= NO_IMPROVEMENT
                    break
                end
            end
        end

        if include_in_population
            current.chromosomes[current.fitness[end][2]] .= best_found[2]
            current.fitness[end] = (best_found[1], current.fitness[end][2])
            # Reorder the chromosomes.
            sort!(current.fitness, rev = sense)
            final_status |= ELITE_IMPROVEMENT
        end
    end

    return final_status
end
