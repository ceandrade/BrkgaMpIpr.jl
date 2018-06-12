################################################################################
# util.jl: utility functions for tests.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jun 11, 2018 by ceandrade
# Last update: Jun 11, 2018 by ceandrade
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

function load_brkga_data(filename::String, brkga_data::BrkgaData)
    tmp = load(filename)

    brkga_data.opt_sense = tmp["opt_sense"]
    brkga_data.chromosome_size = tmp["chromosome_size"]
    brkga_data.population_size = tmp["population_size"]
    brkga_data.elite_size = tmp["elite_size"]
    brkga_data.num_mutants = tmp["num_mutants"]
    brkga_data.num_elite_parents = tmp["num_elite_parents"]
    brkga_data.total_parents = tmp["total_parents"]
    brkga_data.bias = tmp["bias"]
    brkga_data.num_independent_populations = tmp["num_independent_populations"]
    brkga_data.evolutionary_mechanism_on = tmp["evolutionary_mechanism_on"]

    # TODO (ceandrade): list the path relink parameters here.

    # FIXME (ceandrade): the following doesn't work because it tries to
    # load the decoder function from the file. So, we rebuild the instance.
    # brkga_data.problem_instance = tmp["problem_instance"],
    brkga_data.problem_instance =
        TestInstance.Instance(brkga_data.chromosome_size)

    # NOTE (ceandrade): currently, JLD cannot save functions.
    brkga_data.decode! = decode!

    brkga_data.rng = tmp["rng"]
    brkga_data.previous = tmp["previous"]
    brkga_data.current = tmp["current"]
    brkga_data.total_bias_weight = tmp["total_bias_weight"]
    brkga_data.shuffled_individuals = tmp["shuffled_individuals"]
    brkga_data.parents_ordered = tmp["parents_ordered"]
    brkga_data.initialized = tmp["initialized"]
    brkga_data.reset_phase = tmp["reset_phase"]

    # NOTE (ceandrade): currently, JLD cannot save functions.
    if brkga_data.bias == LOGINVERSE
        set_bias_custom_function!(brkga_data, r -> 1.0 / log1p(r))
    elseif brkga_data.bias == LINEAR
        set_bias_custom_function!(brkga_data, r -> 1.0 / r)
    elseif brkga_data.bias == QUADRATIC
        set_bias_custom_function!(brkga_data, r -> r ^ -2.0)
    elseif brkga_data.bias == CUBIC
        set_bias_custom_function!(brkga_data, r -> r ^ -3.0)
    elseif brkga_data.bias == EXPONENTIAL
        set_bias_custom_function!(brkga_data, r -> exp(-r))
    elseif brkga_data.bias == CONSTANT
        set_bias_custom_function!(brkga_data, (::Int64) -> 1.0 / total_parents)
    end

    nothing
end
