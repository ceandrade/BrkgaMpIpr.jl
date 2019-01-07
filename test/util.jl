################################################################################
# util.jl: utility functions for tests.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 11, 2018 by ceandrade
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

using JLD
using FileIO

################################################################################
# Write data file
################################################################################

function write_data(filename::String, data::BrkgaData)
    save(File(format"JLD", filename),
        "opt_sense", brkga_data.opt_sense,
        "chromosome_size", brkga_data.chromosome_size,
        "evolutionary_mechanism_on", brkga_data.evolutionary_mechanism_on,
        "elite_size", brkga_data.elite_size,
        "num_mutants", brkga_data.num_mutants,
        "population_size", brkga_data.params.population_size,
        "num_elite_parents", brkga_data.params.num_elite_parents,
        "total_parents", brkga_data.params.total_parents,
        "bias_type", brkga_data.params.bias_type,
        "num_independent_populations", brkga_data.params.num_independent_populations,
        "pr_number_pairs", brkga_data.params.pr_number_pairs,
        "pr_minimum_distance", brkga_data.params.pr_minimum_distance,
        "pr_type", brkga_data.params.pr_type,
        "pr_selection", brkga_data.params.pr_selection,
        "alpha_block_size", brkga_data.params.alpha_block_size,
        "pr_percentage", brkga_data.params.pr_percentage,
        "problem_instance", brkga_data.problem_instance,
        # NOTE (ceandrade): currently, JLD cannot save functions.
        # decode!::Function
        "rng", brkga_data.rng,
        "previous", brkga_data.previous,
        "current", brkga_data.current,
        # NOTE (ceandrade): currently, JLD cannot save functions.
        # "bias_function", brkga_data.bias_function,
        "total_bias_weight", brkga_data.total_bias_weight,
        "shuffled_individuals", brkga_data.shuffled_individuals,
        "parents_ordered", brkga_data.parents_ordered,
        "initialized", brkga_data.initialized,
        "reset_phase", brkga_data.reset_phase
    )
    nothing
end

################################################################################
# Load data file
################################################################################

function load_brkga_data(filename::String, brkga_data::BrkgaData)
    tmp = load(filename)

    brkga_data.opt_sense = tmp["opt_sense"]
    brkga_data.chromosome_size = tmp["chromosome_size"]
    brkga_data.elite_size = tmp["elite_size"]
    brkga_data.num_mutants = tmp["num_mutants"]
    brkga_data.params.population_size = tmp["population_size"]
    brkga_data.params.num_elite_parents = tmp["num_elite_parents"]
    brkga_data.params.total_parents = tmp["total_parents"]
    brkga_data.params.bias_type = tmp["bias_type"]
    brkga_data.params.num_independent_populations = tmp["num_independent_populations"]
    brkga_data.params.pr_number_pairs = tmp["pr_number_pairs"]
    brkga_data.params.pr_minimum_distance = tmp["pr_minimum_distance"]
    brkga_data.params.pr_type = tmp["pr_type"]
    brkga_data.params.pr_selection = tmp["pr_selection"]
    brkga_data.params.alpha_block_size = tmp["alpha_block_size"]
    brkga_data.params.pr_percentage = tmp["pr_percentage"]

    # FIXME (ceandrade): the following doesn't work because it tries to
    # load the decoder function from the file. So, we rebuild the instance.
    # brkga_data.problem_instance = tmp["problem_instance"],
    brkga_data.problem_instance = Instance(brkga_data.chromosome_size)

    # NOTE (ceandrade): currently, JLD cannot save functions.
    brkga_data.decode! = sum_decode!

    brkga_data.rng = tmp["rng"]
    brkga_data.previous = tmp["previous"]
    brkga_data.current = tmp["current"]
    brkga_data.total_bias_weight = tmp["total_bias_weight"]
    brkga_data.shuffled_individuals = tmp["shuffled_individuals"]
    brkga_data.parents_ordered = tmp["parents_ordered"]
    brkga_data.initialized = tmp["initialized"]
    brkga_data.reset_phase = tmp["reset_phase"]

    # NOTE (ceandrade): currently, JLD cannot save functions.
    if brkga_data.params.bias_type == LOGINVERSE
        set_bias_custom_function!(brkga_data, r -> 1.0 / log1p(r))
    elseif brkga_data.params.bias_type == LINEAR
        set_bias_custom_function!(brkga_data, r -> 1.0 / r)
    elseif brkga_data.params.bias_type == QUADRATIC
        set_bias_custom_function!(brkga_data, r -> r ^ -2.0)
    elseif brkga_data.params.bias_type == CUBIC
        set_bias_custom_function!(brkga_data, r -> r ^ -3.0)
    elseif brkga_data.params.bias_type == EXPONENTIAL
        set_bias_custom_function!(brkga_data, r -> exp(-r))
    elseif brkga_data.params.bias_type == CONSTANT
        set_bias_custom_function!(brkga_data, (::Int64) -> 1.0 / total_parents)
    end

    nothing
end
