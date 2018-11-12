################################################################################
# path_relink.jl: main evolutionary routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 06, 2018 by ceandrade
# Last update: Nov 09, 2018 by ceandrade
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
# Auxiliary structures and functions
################################################################################

"""
    mutable struct Triple

Hold the data structures used to build a candidate chromosome for parallel
decoding on direct path relink.

**THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.**
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

**THIS IS AN INTERNAL DATA STRUCTURE AND IT IS NOT MEANT TO BE USED DIRECTLY.**
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
limited to the `max_end`. **Note:** this function only accept positive numbers,
and all sanity check is disregarded due to performance reasons.

**THIS IS AN INTERNAL FUNCTION AND IT IS NOT MEANT TO BE USED DIRECTLY.**
"""
function find_block_range(block_number::Int64, block_size::Int64,
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

**Note:** this function only accept positive numbers, and all sanity and
bounds check is disregarded due to performance reasons.

**THIS IS AN INTERNAL FUNCTION AND IT IS NOT MEANT TO BE USED DIRECTLY.**
"""
@inline function swap!(x::Array{T, 1}, pos1::Int64, pos2::Int64) where T
    @inbounds x[pos1], x[pos2] = x[pos2], x[pos1]
    nothing
end

################################################################################
# Direct path relinking method
################################################################################

"""
    function direct_path_relink!(brkga_data::BrkgaData,
                                 population_index::Int64,
                                 chr1_index::Int64,
                                 chr2_index::Int64,
                                 affect_solution::Function,
                                 block_size::Int64,
                                 max_time::Int64,
                                 percentage::Float64
        )::Tuple{Float64, Array{Float64, 1}}

Perform the direct path relinking, changing each allele or block
of alleles of base chromosome for the correspondent one in the guide
chromosome.

The API will call `decode!()` function, in `BrkgaData`, always with
`writeback = false`. The reason is that if the decoder rewrites the chromosome,
the path between solutions is lost and inadvertent results may come up.
Note that at the end of the path relinking, the method calls the decoder with
`writeback = true` in the best chromosome found to guarantee that this
chromosome is re-written to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).

**THIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS
CALLED FROM THE `path_relink()` FUNCTION.** Due to this reason, this method
**DOES NOT** perform health checks on the arguments.

# Arguments
- `brkga_data::BrkgaData`: the BRKGA data.

- `population_index::Int64`: the population from where the chromosomes will be
  analized.

- `chr1_index::Int64` and `chr2_index::Int64`: two valid indices from
  chromosomes of the given population.

- `affect_solution::Function`: function that takes two partial chromosomes /
  block of genes `BLOCK1` and `BLOCK2` and checks whether changing the keys from
  `BLOCK1` to `BLOCK2` affects the solution. For instance, suppose that the
  alleles/keys are used as threshold such that values > 0.5 activate a feature.
  Suppose we have `BLOCK1 = [0.3, 0.4, 0.1]` and `BLOCK2 = [0.4, 0.1, 0.2]`.
  Since all values are below 0.5, changing the keys from `BLOCK1` to `BLOCK2`
  does not chage the solution, and therefore, we can drop such change (and
  subsequentely decoding). The blocks can hold only one key/allele, sequential
  key blocks, of even the whole chromosome. `affect_solution` takes two
  views/subarrays. The function **must have** the following signature

        `affect_solution(block_1::SubArray{Float64, 1},
                         block_2::SubArray{Float64, 1})::Float64`

  **Note: this function depends on the problem structure and how the
  keys/alleles are used.**

- `block_size::Int64`: (posite) number of alleles to be exchanged at once in
  each iteration. If `block_size == 1`, the traditional path relinking is
  performed.

- `max_time::Int64`: abort path-relinking when reach `max_time`.
       If `max_time <= 0`, no limit is imposed. Given in seconds.

- `percentage::Float64`: define the size, in percentage, of the path to
       build. Range [0, 1].

# Returns

- `Tuple{Float64, Array{Float64, 1}}`: the best pair (fitness,   chromosome)
  found during the relinking. If the relink is not possible due to homogeneity,
  `-Inf` returns in case of maximization, and `Inf` in case of minimization.

"""
function direct_path_relink!(brkga_data::BrkgaData,
                             population_index::Int64,
                             chr1_index::Int64,
                             chr2_index::Int64,
                             affect_solution::Function,
                             block_size::Int64,
                             max_time::Int64,
                             percentage::Float64
    )::Tuple{Float64, Array{Float64, 1}}

    PR_START_TIME = time()

    bd = brkga_data

    best_chr_found = Array{Float64, 1}(undef, bd.chromosome_size)
    best_fitness_found::Float64 = (bd.opt_sense == MAXIMIZE) ? -Inf : Inf

    NUM_BLOCKS = cld(bd.chromosome_size, block_size)
    PATH_SIZE = Int64(floor(percentage * NUM_BLOCKS))
    remaining_blocks = collect(1:NUM_BLOCKS)

    # Allocate memory for the candidates.
    candidates_base = Array{Triple, 1}(undef, NUM_BLOCKS)
    candidates_guide = Array{Triple, 1}(undef, NUM_BLOCKS)

    # References for the base and guide chromosomes.
    base = bd.current[population_index].chromosomes[chr1_index]
    guide = bd.current[population_index].chromosomes[chr2_index]

    # **NOTE:** two loops are faster than one according to `@benchmark`,
    # probably because of cache racing.
    @inbounds Threads.@threads for i in 1:NUM_BLOCKS
        candidates_base[i] = Triple(copy(base))
    end
    @inbounds Threads.@threads for i in 1:NUM_BLOCKS
        candidates_guide[i] = Triple(copy(guide))
    end

    # Holds the original keys.
    old_keys = Array{Float64, 1}(undef, brkga_data.chromosome_size)

    SENSE = (bd.opt_sense == MAXIMIZE)
    iterations = 1
    while !isempty(remaining_blocks)
        state = iterate(remaining_blocks)

        i = 1
        while state !== nothing && !isempty(remaining_blocks)
            it_block_idx = state[2] - 1
            BLOCK_RANGE = find_block_range(remaining_blocks[it_block_idx],
                                           block_size, length(guide))
            BLOCK1 = view(candidates_base[i].chr, BLOCK_RANGE)
            BLOCK2 = view(guide, BLOCK_RANGE)

            if !affect_solution(BLOCK1, BLOCK2)
                deleteat!(remaining_blocks, it_block_idx)
                continue
            end

            # Save the former keys before...
            old_keys[BLOCK_RANGE] = candidates_base[i].chr[BLOCK_RANGE]

            # ... copy the keys from the guide solution.
            candidates_base[i].chr[BLOCK_RANGE] = guide[BLOCK_RANGE]

            candidates_base[i].block_index = it_block_idx
            state = iterate(remaining_blocks, state[2])
            i += 1
        end

        if i == 1 && isempty(remaining_blocks)
            break
        end

        Threads.@threads for i in 1:length(remaining_blocks)
            candidates_base[i].fitness = (SENSE) ? -Inf : Inf
            if time() - PR_START_TIME > max_time
                continue
            end

            candidates_base[i].fitness =
                bd.decode!(candidates_base[i].chr, bd.problem_instance, false)
        end

        # Locate the best candidate.
        best_value::Float64 = (SENSE) ? -Inf : Inf
        best_index = 1
        best_block_index = 1

        @inbounds for i in 1:length(remaining_blocks)
            if ((best_value < candidates_base[i].fitness && SENSE) ||
                (best_value > candidates_base[i].fitness && !SENSE))
                best_value = candidates_base[i].fitness
                best_index = i
                best_block_index = candidates_base[i].block_index
            end
        end

        # Hold it, if it is the best found until now.
        if ((SENSE && best_fitness_found < candidates_base[best_index].fitness)
            ||
            (!SENSE && best_fitness_found > candidates_base[best_index].fitness))
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
        if iterations == PATH_SIZE || time() - PR_START_TIME > max_time
            break
        end
    end

    if best_fitness_found == Inf || best_fitness_found == -Inf
        best_chr_found = Array{Float64, 1}()
    end
    return (best_fitness_found, best_chr_found)
end

################################################################################
# Permutation-based path relinking method
################################################################################

"""
    function permutation_based_path_relink!(brkga_data::BrkgaData,
                                            population_index::Int64,
                                            chr1_index::Int64,
                                            chr2_index::Int64,
                                            affect_solution::Function,
                                            block_size::Int64,
                                            max_time::Int64,
                                            percentage::Float64
        )::Tuple{Float64, Array{Float64, 1}}

Performs the permutation-based path relinking. In this method, the permutation
induced by the keys in the guide solution is used to change the order of the
keys in the permutation induced by the base solution.

The API will call `decode!()` function, in `BrkgaData`, always with
`writeback = false`. The reason is that if the decoder rewrites the chromosome,
the path between solutions is lost and inadvertent results may come up.
Note that at the end of the path relinking, the method calls the decoder with
`writeback = true` in the best chromosome found to guarantee that this
chromosome is re-written to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).

**THIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS
CALLED FROM THE `path_relink()` FUNCTION.** Due to this reason, this method
**DOES NOT** perform health checks on the arguments.

# Arguments
- `brkga_data::BrkgaData`: the BRKGA data.

- `population_index::Int64`: the population from where the chromosomes will be
  analized.

- `chr1_index::Int64` and `chr2_index::Int64`: two valid indices from
  chromosomes of the given population.

- `affect_solution::Function`: not used in this function but kept to API
  compatibility.

- `block_size::Int64`: not used in this function but kept to API compatibility.

- `max_time::Int64`: abort path-relinking when reach `max_time`.
       If `max_time <= 0`, no limit is imposed. Given in seconds.

- `percentage::Float64`: define the size, in percentage, of the path to
       build. Range [0, 1].

# Returns

- `Tuple{Float64, Array{Float64, 1}}`: the best pair (fitness,   chromosome)
  found during the relinking. If the relink is not possible due to homogeneity,
  `-Inf` returns in case of maximization, and `Inf` in case of minimization.

"""
function permutation_based_path_relink!(brkga_data::BrkgaData,
                                        population_index::Int64,
                                        chr1_index::Int64,
                                        chr2_index::Int64,
                                        affect_solution::Function,
                                        block_size::Int64,
                                        max_time::Int64,
                                        percentage::Float64
    )::Tuple{Float64, Array{Float64, 1}}

     PR_START_TIME = time()

    bd = brkga_data

    best_chr_found = Array{Float64, 1}(undef, bd.chromosome_size)
    best_fitness_found::Float64 = (bd.opt_sense == MAXIMIZE) ? -Inf : Inf

    PATH_SIZE = Int64(floor(percentage * bd.chromosome_size))
    remaining_indices = collect(1:bd.chromosome_size)

    # Allocate memory for the candidates.
    candidates_base = Array{DecodeStruct, 1}(undef, bd.chromosome_size)
    candidates_guide = Array{DecodeStruct, 1}(undef, bd.chromosome_size)

    # References for the base and guide chromosomes.
    base = bd.current[population_index].chromosomes[chr1_index]
    guide = bd.current[population_index].chromosomes[chr2_index]

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

    SENSE = (bd.opt_sense == MAXIMIZE)
    iterations = 1
    while !isempty(remaining_indices)
        # it_idx = start(remaining_indices)
        state = iterate(remaining_indices)

        i = 1
        # while !done(remaining_indices, it_idx)
        while state !== nothing && !isempty(remaining_indices)
            it_idx = state[2] - 1
            position_in_base = base_indices[remaining_indices[it_idx]]
            position_in_guide = guide_indices[remaining_indices[it_idx]]

            if position_in_base == position_in_guide
                deleteat!(remaining_indices, it_idx)
                continue
            end

            candidates_base[i].key_index = it_idx
            candidates_base[i].pos1 = position_in_base
            candidates_base[i].pos2 = position_in_guide
            candidates_base[i].fitness = (SENSE) ? -Inf : Inf

            # _, it_idx = next(remaining_indices, it_idx)
            state = iterate(remaining_indices, state[2])
            i += 1
        end

        if i == 1 && isempty(remaining_indices)
            break
        end

        # Decode the candidates.
        Threads.@threads for i in 1:length(remaining_indices)
            if time() - PR_START_TIME > max_time
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
        best_value::Float64 = (SENSE) ? -Inf : Inf
        best_index = 1
        best_key_index = 1

        @inbounds for i in 1:length(remaining_indices)
            if ((best_value < candidates_base[i].fitness && SENSE) ||
                (best_value > candidates_base[i].fitness && !SENSE))
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
        if ((SENSE && best_fitness_found < candidates_base[best_index].fitness)
            ||
            (!SENSE && best_fitness_found > candidates_base[best_index].fitness))
            best_fitness_found = candidates_base[best_index].fitness
            best_chr_found[:] = candidates_base[best_index].chr
        end

        # Swap roles.
        base, guide = guide, base
        candidates_base, candidates_guide = candidates_guide, candidates_base
        deleteat!(remaining_indices, best_key_index)

        iterations += 1
        if iterations == PATH_SIZE || time() - PR_START_TIME > max_time
            break
        end
    end

    if best_fitness_found == Inf || best_fitness_found == -Inf
        best_chr_found = Array{Float64, 1}()
    end
    return (best_fitness_found, best_chr_found)
end

################################################################################
# Main path relinking method
################################################################################

"""
    function permutation_based_path_relink!(brkga_data::BrkgaData,
                                            population_index::Int64,
                                            chr1_index::Int64,
                                            chr2_index::Int64,
                                            affect_solution::Function,
                                            block_size::Int64,
                                            max_time::Int64,
                                            percentage::Float64
        )::Tuple{Float64, Array{Float64, 1}}

Performs the permutation-based path relinking. In this method, the permutation
induced by the keys in the guide solution is used to change the order of the
keys in the permutation induced by the base solution.

The API will call `decode!()` function, in `BrkgaData`, always with
`writeback = false`. The reason is that if the decoder rewrites the chromosome,
the path between solutions is lost and inadvertent results may come up.
Note that at the end of the path relinking, the method calls the decoder with
`writeback = true` in the best chromosome found to guarantee that this
chromosome is re-written to reflect the best solution found.

This method is a multi-thread implementation. Instead of to build and
decode each chromosome one at a time, the method builds a list of
candidates, altering the alleles/keys according to the guide solution,
and then decode all candidates in parallel. Note that
`O(chromosome_size^2 / block_size)` additional memory is necessary to build
the candidates, which can be costly if the `chromosome_size` is very large.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).

**THIS IS AN INTERNAL METHOD AND IT IS NOT MEANT TO BE USED DIRECTLY. IT IS
CALLED FROM THE `path_relink()` FUNCTION.** Due to this reason, this method
**DOES NOT** perform health checks on the arguments.

# Arguments
- `brkga_data::BrkgaData`: the BRKGA data.

- `population_index::Int64`: the population from where the chromosomes will be
  analized.

- `chr1_index::Int64` and `chr2_index::Int64`: two valid indices from
  chromosomes of the given population.

- `affect_solution::Function`: not used in this function but kept to API
  compatibility.

- `block_size::Int64`: not used in this function but kept to API compatibility.

- `max_time::Int64`: abort path-relinking when reach `max_time`.
       If `max_time <= 0`, no limit is imposed. Given in seconds.

- `percentage::Float64`: define the size, in percentage, of the path to
       build. Range [0, 1].

# Returns

- `Tuple{Float64, Array{Float64, 1}}`: the best pair (fitness,   chromosome)
  found during the relinking. If the relink is not possible due to homogeneity,
  `-Inf` returns in case of maximization, and `Inf` in case of minimization.

"""
function path_relink!(brkga_data::BrkgaData,
                      compute_distance::Function,
                      affect_solution::Function,
                      minimum_distance::Float64,
                      pr_type::PathRelinkingType,
                      pr_selection::PathRelinkingSelection,
                      block_size::Int64,
                      max_time::Int64,
                      percentage::Float64
    )::Bool

    return true
end
