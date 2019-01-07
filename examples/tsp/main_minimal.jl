################################################################################
# main_minimal.jl: minimal script for call BRKGA algorithms to solve
#                  instances of the Traveling Salesman Problem.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Dec 28, 2018 by ceandrade
# Last update: Jan 07, 2019 by ceandrade
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

import Random
using BrkgaMpIpr

include("tsp_instance.jl")
include("tsp_decoder.jl")

########################################
# Read arguments and Instance
########################################

if length(ARGS) < 4
    println("Usage: julia main_minimal.jl <seed> <config-file> " *
            "<num-generations> <tsp-instance-file>")
    exit(1)
end

seed = parse(Int64, ARGS[1])
configuration_file = ARGS[2]
num_generations = parse(Int64, ARGS[3])
instance_file = ARGS[4]

println("Reading data...")
instance = TSP_Instance(instance_file)

########################################
# Build the BRKGA data structures and initialize
########################################

println("Building BRKGA data and initializing...")
(brkga_data, control_params) = build_brkga(
    instance, tsp_decode!, MINIMIZE, seed, instance.num_nodes,
    configuration_file
)

# NOTE: don't forget to initialize the algorithm.
initialize!(brkga_data)

########################################
# Find good solutions / evolve
########################################

println("Evolving $num_generations generations...")
evolve!(brkga_data, num_generations)

best_solution = get_best_fitness(brkga_data)
@show best_solution
