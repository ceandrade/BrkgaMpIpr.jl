################################################################################
# support.jl: initialization, reset, and individual migration routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 26, 2018 by ceandrade
# Last update: Mar 27, 2018 by ceandrade
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

"""
    initialize!(brkga_data::BrkgaData)

Initialize the populations and others data structures of the BRKGA. If an
initial population is supplied, this method completes the remaining individuals,
if they do not exist.

**THIS METHOD MUST BE CALLED BEFORE ANY OPTIMIZATIOM METHODS.**

This method also performs the initial decoding of the chromosomes. Therefore,
depending on the decoder implementation, this can take a while, and the user may
want to time such procedure in his/her experiments.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).
"""
function initialize!(brkga_data::BrkgaData)
    bd = brkga_data  # Just an short alias.

    # If we have warmstaters, complete the population if necessary.
    # Note that it is done only in the true initialization.
    pop_start = 1
    if length(bd.current) > 0 && !bd.reset_phase
        population = bd.current[1]
        for i = (length(population.chromosomes) + 1):bd.population_size
            push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
        end
        population.fitness = Array{Tuple{Float64, Int64}, 1}(bd.population_size)
        pop_start = 2

    elseif length(bd.current) == 0
        bd.current = Array{Population, 1}(bd.num_independent_populations)
        bd.previous = Array{Population, 1}(bd.num_independent_populations)
    end

    # Build the remaining populations and associated data structures.
    for i = pop_start:bd.num_independent_populations
        # If no reset, allocate memory.
        if !bd.reset_phase
            population = Population()
            for j = 1:bd.population_size
                push!(population.chromosomes, rand(bd.rng, bd.chromosome_size))
            end
            population.fitness =
                Array{Tuple{Float64, Int64}, 1}(bd.population_size)
            brkga_data.current[i] = population

        # For reset phase, only generates new keys. No new memory allocation
        # and much faster than copy/deepcopy.
        else
            for chr in brkga_data.current[i].chromosomes
                for j = 1:length(chr)
                    chr[j] = rand(bd.rng)
                end
            end
        end
    end

    # Perform initial decoding. It may take a while.
    for population in bd.current
        Threads.@threads for i in eachindex(population.chromosomes)
            value = bd.decode!(population.chromosomes[i], bd.problem_instance)
            population.fitness[i] = (value, i)
        end
        sort!(population.fitness, rev = bd.opt_sense == MAXIMIZE)
    end

    # Copy the data to previous populations.
    # NOTE (ceandrade) During reset phase, copying item by item maybe faster
    # than deepcoping (which allocates new memory).
    bd.previous = deepcopy(bd.current)

    bd.initialized = true;
    bd.reset_phase = false;
    nothing
end

################################################################################

"""
    reset!(brkga_data::BrkgaData)

Reset all populations with brand new keys. All warm start solutions provided
by `set_initial_population!()` are discarded.

**NOTE:** as it is in `evolve()`, the decoding is done in parallel using
threads, and the user **must guarantee that the decoder is THREAD-SAFE.**
If such property cannot be held, we suggest using single thread by setting the
environmental variable `JULIA_NUM_THREADS = 1`
(see https://docs.julialang.org/en/stable/manual/parallel-computing).
"""
function reset!(brkga_data::BrkgaData)
    if !brkga_data.initialized
        error("the algorithm hasn't been initialized. Call initialize!() before reset!()")
    end
    brkga_data.reset_phase = true
    initialize!(brkga_data)
end
