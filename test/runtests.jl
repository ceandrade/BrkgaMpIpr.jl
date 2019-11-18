#!/usr/bin/env julia
################################################################################
# runtests.jl: unit tests for BrkgaMpIpr
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
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

using BrkgaMpIpr
using JLD
using Test
using Random

################################################################################
# General fixtures for testing
################################################################################

include("test_instance.jl")
include("test_decoders.jl")

const CHROMOSOME_SIZE = 100
instance = Instance(CHROMOSOME_SIZE)

# Makes easy to change specific position on the parameters vector below.
const param_names = ["instance", "decode!", "opt_sense", "seed", "chr_size",
                     "brkga_params", "evolutionary_mechanism_on"]

# Reverse index.
const param_index = Dict([v => i for (i, v) in enumerate(param_names)])

# Holds the parameters to build new BrkgaData.
const default_param_values = Array{Any, 1}(undef, length(param_names))

default_brkga_params = BrkgaParams()
default_brkga_params.population_size = 10
default_brkga_params.elite_percentage = 0.3
default_brkga_params.mutants_percentage = 0.1
default_brkga_params.num_elite_parents = 1
default_brkga_params.total_parents = 2
default_brkga_params.bias_type = LOGINVERSE
default_brkga_params.num_independent_populations = 3
default_brkga_params.pr_number_pairs = 0
default_brkga_params.pr_minimum_distance = 0.0
default_brkga_params.pr_type = DIRECT
default_brkga_params.pr_selection = BESTSOLUTION
default_brkga_params.alpha_block_size = 1.0
default_brkga_params.pr_percentage = 1.0

# Some default parameters.
default_param_values[param_index["instance"]] = instance
default_param_values[param_index["decode!"]] = sum_decode!
default_param_values[param_index["opt_sense"]] = MAXIMIZE
default_param_values[param_index["seed"]] = 2700001
default_param_values[param_index["chr_size"]] = CHROMOSOME_SIZE
default_param_values[param_index["brkga_params"]] = default_brkga_params
default_param_values[param_index["evolutionary_mechanism_on"]] = true

################################################################################
# Functions
################################################################################

@time begin
    print(">> Testing types operations...\n")
    include("types_tests.jl")

    print(">> Testing types I/O operations...\n")
    include("types_io_tests.jl")

    print("\n>> Testing BRKGA data building and initialization...\n")
    include("building_tests.jl")

    print("\n>> Testing support methods...\n")
    include("support_tests.jl")

    print("\n>> Testing population manipulation methods...\n")
    include("population_manipulation_tests.jl")

    print("\n>> Testing distance functions...\n")
    include("distance_functions_tests.jl")

    print("\n>> Testing evolutionary methods (it may take a while)...\n")
    include("evolution_tests.jl")

    print("\n>> Testing path relink methods (it may take a while)...\n")
    # include("path_relink_tests.jl")
end
