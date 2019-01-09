################################################################################
# main_complete.jl: comprehensive script for BRKGA-MP-IPR experiments
#                   using Julia.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jan 07, 2019 by ceandrade
# Last update: Jan 09, 2019 by ceandrade
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

push!(LOAD_PATH, "../..")
push!(LOAD_PATH, ".")

using BrkgaMpIpr

import Base: parse
import Dates
import DocOpt
import Random
using Printf

include("tsp_instance.jl")
include("greedy_tour.jl")
include("tsp_decoder.jl")

################################################################################
# Enumerations and constants
################################################################################

"""
    @enum StopRule

Controls stop criteria. Stops either when:
- a given number of `GENERATIONS` is given;
- or a `TARGET` value is found;
- or no `IMPROVEMENT` is found in a given number of iterations.
"""
@enum StopRule begin
    GENERATIONS = 0
    TARGET = 1
    IMPROVEMENT = 2
end

"""
    function parse(::Type{StopRule}, value::String)::StopRule

Parse `value` into a `StopRule`.
"""
function parse(::Type{StopRule}, value::String)::StopRule
    local_value = uppercase(strip(value)[1])
    if local_value == 'G'
        return GENERATIONS
    elseif local_value == 'T'
        return TARGET
    elseif local_value == 'I'
        return IMPROVEMENT
    end

    throw(ArgumentError("cannot parse $value as StopRule"))
end

################################################################################
# Main function
################################################################################

"""
    function main(args)

Proceed with the optimization. Create to avoid spread `global` keywords around
the code.
"""
function main(args)
    configuration_file = args["<config_file>"]
    instance_file = args["<instance_file>"]
    seed = parse(Int64, args["<seed>"])
    stop_rule = parse(StopRule, args["<stop_rule>"])

    if stop_rule == TARGET
        stop_argument = parse(Float64, args["<stop_arg>"])
    else
        stop_argument = parse(Int64, args["<stop_arg>"])
    end

    maximum_time = parse(Float64, args["<max_time>"])
    if maximum_time <= 0.0
        error("Maximum time must be larger than 0.0. Given $maximum_time.")
    end

    perform_evolution = !args["--no_evolution"]

    ########################################
    # Load config file and show basic info.
    ########################################

    brkga_params, control_params = load_configuration(configuration_file)

    print("""
    ------------------------------------------------------
    > Experiment started at $(Dates.now())
    > Instance: $instance_file
    > Configuration: $configuration_file
    > Algorithm Parameters:
    """)

    if !perform_evolution
        println(">    - Simple multi-start: on (no evolutionary operators)")
    else
        output_string = ""
        for field in fieldnames(BrkgaParams)
            output_string *= ">  - $field $(getfield(brkga_params, field))\n"
        end
        for field in fieldnames(ExternalControlParams)
            output_string *= ">  - $field $(getfield(control_params, field))\n"
        end
        print(output_string)
        println("""
        > Seed: $seed
        > Stop rule: $stop_rule
        > Stop argument: $stop_argument
        > Maximum time (s): $maximum_time
        > Number of parallel threads for decoding: $(Threads.nthreads())
        ------------------------------------------------------""")
    end

    ########################################
    # Load instance and adjust BRKGA parameters
    ########################################

    println("\n[$(Dates.Time(Dates.now()))] Reading TSP data...")

    instance = TSP_Instance(instance_file)
    println("Number of nodes: $(instance.num_nodes)")

    println("\n[$(Dates.Time(Dates.now()))] Generating initial tour...")

    # Generate a greedy solution to be used as warm start for BRKGA.
    initial_cost, initial_tour = greedy_tour(instance)
    println("Initial cost: $initial_cost")

    ########################################
    # Build the BRKGA data structures and initialize
    ########################################

    println("\n[$(Dates.Time(Dates.now()))] Building BRKGA data...")

    # Usually, it is a good idea to set the population size
    # proportional to the instance size.
    brkga_params.population_size = min(brkga_params.population_size,
                                       10 * instance.num_nodes)
    println("New population size: $(brkga_params.population_size)")

    # Chromosome size is the number of nodes.
    # Each chromosome represents a permutation of nodes.
    brkga_data = build_brkga(instance, tsp_decode!, MINIMIZE, seed,
                            instance.num_nodes, brkga_params, perform_evolution)

    # To inject the initial tour, we need to create chromosome representing that
    # solution. First, we create a set of keys to be used in the chromosome.
    Random.seed!(seed)
    keys = sort(rand(instance.num_nodes))

    # Then, we visit each node in the tour and assign to it a key.
    initial_chromosome = zeros(instance.num_nodes)
    for i in 1:instance.num_nodes
        initial_chromosome[initial_tour[i]] = keys[i]
    end

    # Inject the warm start solution in the initial population.
    set_initial_population!(brkga_data, [initial_chromosome])

    # NOTE: don't forget to initialize the algorithm.
    println("\n[$(Dates.Time(Dates.now()))] Initializing BRKGA data...")
    initialize!(brkga_data)

    ########################################
    # Warm start the script/code
    ########################################

    # To make sure we are timing the runs correctly, we run some warmup
    # iterations with bogus data. Warmup is always recommended for script
    # languages. Here, we call the most used methods.
    println("\n[$(Dates.Time(Dates.now()))] Warming up...")

    bogus_data = deepcopy(brkga_data)
    evolve!(bogus_data, 2)
    path_relink!(bogus_data, (x, y) -> 1.0, (x, y) -> true, 0, 0.5,
                 brkga_params.pr_type, brkga_params.pr_selection, 1, 10.0, 1.0)
    get_best_fitness(brkga_data)
    get_best_chromosome(brkga_data)

    bogus_data = nothing

    ########################################
    # Evolving
    ########################################

    println("\n[$(Dates.Time(Dates.now()))] Evolving...")
    println("* Iteration | Cost | CurrentTime")

    best_cost = Inf
    best_chromosome = initial_chromosome

    iteration  = 0
    last_update_time = 0.0
    last_update_iteration = 0
    large_offset = 0
    path_relink_time = 0.0
    num_path_relink_calls = 0
    num_homogenities = 0
    num_best_improvements = 0
    num_elite_improvements = 0
    run = true
    start_time = time()

    # Main optimization loop. We evolve one generation at time,
    # keeping track of all changes during such process.
    while run
        iteration += 1

        # Evolves one iteration.
        evolve!(brkga_data)

        # Checks the current results and holds the best.
        fitness = get_best_fitness(brkga_data)
        if fitness < best_cost
            last_update_time = time() - start_time
            update_offset = iteration - last_update_iteration

            if large_offset < update_offset
                large_offset = update_offset
            end

            last_update_iteration = iteration
            best_cost = fitness
            best_chromosome = get_best_chromosome(brkga_data)

            @printf("* %d | %.0f | %.2f \n", iteration, best_cost,
                    last_update_time)
        end

        iter_without_improvement = iteration - last_update_iteration

        # Here, we call the path relink when the algorithm gets stuck for
        # `exchange_interval` iterations. Obvisualy, we can use many other ways
        # of hybridization.
        if control_params.exchange_interval > 0 &&
           iter_without_improvement > 0 &&
           (iter_without_improvement % control_params.exchange_interval == 0)

            println("Performing path relink at $iteration...")
            num_path_relink_calls += 1

            pr_now = time()
            result = path_relink!(
                brkga_data, kendall_tau_distance, affect_solution_kendall_tau,
                brkga_params.pr_number_pairs,
                brkga_params.pr_minimum_distance,
                brkga_params.pr_type,
                brkga_params.pr_selection,
                1, #block_size doesn't not matter for permutation.
                maximum_time - (time() - start_time),
                brkga_params.pr_percentage
            )

            pr_time = time() - pr_now
            path_relink_time += pr_time

            if result == TOO_HOMOGENEOUS
                num_homogenities += 1
                println("- Populations are too too homogeneous | " *
                        "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == NO_IMPROVEMENT
                println("- No improvement found | " *
                        "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == ELITE_IMPROVEMENT
                num_elite_improvements += 1
                println("- Improvement on the elite set but " *
                        "not in the best individual | " *
                        "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == BEST_IMPROVEMENT
                num_best_improvements += 1
                fitness = get_best_fitness(brkga_data)
                println("- Best individual improvement: $fitness | " *
                        "Elapsed time: $(@sprintf("%.2f", pr_time))")
                if fitness < best_cost
                    last_update_time = time() - start_time
                    update_offset = iteration - last_update_iteration

                    if large_offset < update_offset
                        large_offset = update_offset
                    end

                    last_update_iteration = iteration
                    best_cost = fitness
                    best_chromosome = get_best_chromosome(brkga_data)

                    @printf("* %d | %.0f | %.2f \n", iteration, best_cost,
                            last_update_time)
                end
            end
        end

        # Check stop criteria.
        run = !(
            (time() - start_time > maximum_time) ||
            (stop_rule == GENERATIONS && Float64(iteration) == stop_argument) ||
            (stop_rule == IMPROVEMENT &&
             Float64(iter_without_improvement) >= stop_argument) ||
            (stop_rule == TARGET && best_cost <= stop_argument)
        )
    end
    total_elapsed_time = time() - start_time
    total_num_iterations = iteration

    println("[$(Dates.Time(Dates.now()))] End of optimization")
    print("\nTotal number of iterations: $total_num_iterations")
    print("\nLast update iteration: $last_update_iteration")
    @printf("\nTotal optimization time: %.2f", total_elapsed_time)
    @printf("\nLast update time: %.2f", last_update_time)
    print("\nLarge number of iterations between improvements: $large_offset")

    @printf("\nTotal path relink time: %.2f", path_relink_time)
    print("\nTotal path relink calls: $num_path_relink_calls")
    print("\nNumber of homogenities: $num_homogenities")
    print("\nImprovements in the elite set: $num_elite_improvements")
    print("\nBest individual improvements: $num_best_improvements")

    tour = Array{Tuple{Float64, Int64}}(undef, instance.num_nodes)
    for (index, key) in enumerate(best_chromosome)
        tour[index] = (key, index)
    end
    sort!(tour)

    print("\n\n% Best tour cost: $(@sprintf("%.0f", best_cost))")
    print("\n% Best tour: ")
    for (_, node) in tour
        print("$node ")
    end

    println("\n\nInstance,Seed,NumNodes,TotalIterations,TotalTime," *
            "TotalPRTime,PRCalls,NumHomogenities,NumPRImprovElite," *
            "NumPrImprovBest,LargeOffset,LastUpdateIteration,LastUpdateTime," *
            "Cost")
    print("$(basename(instance_file))," *
           "$seed,$(instance.num_nodes),$total_num_iterations," *
          "$(@sprintf("%.2f", total_elapsed_time))," *
          "$(@sprintf("%.2f", path_relink_time))," *
          "$num_path_relink_calls,$num_homogenities,$num_elite_improvements," *
          "$num_best_improvements,$large_offset,$last_update_iteration," *
          "$(@sprintf("%.2f", last_update_time))," *
          "$(@sprintf("%.0f", best_cost))")
    nothing
end

################################################################################
# Parse and validate command-line arguments, and call main function
################################################################################

doc = """
Usage:
  main_complete.jl -c <config_file> -s <seed> -r <stop_rule> -a <stop_arg> -t <max_time> -i <instance_file> [--no_evolution]
  main_complete.jl (-h | --help)

Options:

    -c <config_file>    Text file with the BRKGA-MP-IPR parameters.

    -s <seed>           Seed for the random number generator.

    -r <stop_rule>      Stop rule where:
                        - (G)enerations: number of evolutionary generations.
                        - (I)terations: maximum number of generations without
                          improvement in the solutions.
                        - (T)arget: runs until obtains the target value.

    -a <stop_arg>       Argument value for '-r'.

    -t <max_time>       Maximum time in seconds.

    -i <instance_file>  Instance file.

    --no_evolution      If supplied, no evolutionary operators are applied. So,
                        the algorithm becomes a simple multi-start algorithm.

    -h --help           Produce help message.
"""
args = DocOpt.docopt(doc)
main(args)
