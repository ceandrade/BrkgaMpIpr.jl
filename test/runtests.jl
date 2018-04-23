#!/usr/bin/env julia
################################################################################
# runtests.jl: unit tests for BrkgaMpIpr
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Apr 23, 2018 by ceandrade
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
using Base.Test

################################################################################
# General fixtures for testing
################################################################################

include("TestInstance.jl")
include("TestDecoder.jl")

using TestInstance
using TestDecoder

const CHROMOSOME_SIZE = 100

instance = Instance(CHROMOSOME_SIZE)

# Makes easy to change espeficif position on the parameters vector below.
const param_names = ["instance", "decode!", "opt_sense", "seed", "chr_size",
                     "pop_size", "elite_percentage", "mutants_percentage",
                     "evolutionary_mechanism_on", "num_elite_parents",
                     "total_parents", "bias", "num_independent_populations"]

# Reverse index.
const param_index = Dict([v => i for (i, v) in enumerate(param_names)])

# Holds the parameters to build new BrkgaData.
param_values = Array{Any, 1}(length(param_names))

# Some default parameters.
param_values[param_index["instance"]] = instance
param_values[param_index["decode!"]] = decode!
param_values[param_index["opt_sense"]] = MAXIMIZE
param_values[param_index["seed"]] = 2700001
param_values[param_index["chr_size"]] = CHROMOSOME_SIZE
param_values[param_index["pop_size"]] = 10
param_values[param_index["elite_percentage"]] = 0.3
param_values[param_index["mutants_percentage"]] = 0.1
param_values[param_index["evolutionary_mechanism_on"]] = true
param_values[param_index["num_elite_parents"]] = 1
param_values[param_index["total_parents"]] = 2
param_values[param_index["bias"]] = LOGINVERSE
param_values[param_index["num_independent_populations"]] = 3

# Used to restore original param_values.
const default_param_values = copy(param_values)

################################################################################
# Functions
################################################################################

print(">> Testing types and their I/O operations...\n")
@time include("types_io_tests.jl")

print("\n>> Testing BRKGA data building and initialization...\n")
@time include("building_tests.jl")

print("\n>> Testing support methods...\n")
@time include("support_tests.jl")

print("\n>> Testing evolutionary methods (it may take a while)...\n")
@time include("evolution_tests.jl")
