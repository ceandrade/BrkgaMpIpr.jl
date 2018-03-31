################################################################################
# types_io.jl: Input/output/parsing methods for internal data strucutures.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 24, 2018 by ceandrade
# Last update: Mar 30, 2018 by ceandrade
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

import Base: parse

################################################################################

"""
    parse(::Type{BiasFunction}, value::String)::BiasFunction

Parse `value` and return a valid `BiasFunction` enumeration.

# Throws
- `ArgumentError`: in case the bias description does not match.
"""
function parse(::Type{BiasFunction}, value::String)::BiasFunction
    value = uppercase(strip(value))
    if value == "CONSTANT"
        return CONSTANT
    elseif value == "CUBIC"
        return CUBIC
    elseif value == "EXPONENTIAL"
        return EXPONENTIAL
    elseif value == "LINEAR"
        return LINEAR
    elseif value == "LOGINVERSE"
        return LOGINVERSE
    elseif value == "QUADRATIC"
        return QUADRATIC
    end

    throw(ArgumentError("cannot parse $value as BiasFunction"))
end

################################################################################

"""
"""
function write_configuration(brkga_data::BrkgaData,
                             external_params::ExternalControlParams,
                             filename::String)

    elite_percentage = brkga_data.elite_size / brkga_data.population_size;
    mutants_percentage = brkga_data.num_mutants / brkga_data.population_size;

    # TODO (ceandrade): implement the path relink parameters.
    open(filename, "w") do file
        write(file,
              "population_size $(brkga_data.population_size)\n",
              "elite_percentage $(elite_percentage)\n",
              "mutants_percentage $(mutants_percentage)\n",
              "mutants_percentage $(mutants_percentage)\n",
              "elite_parents $(brkga_data.num_elite_parents)\n",
              "total_parents $(brkga_data.total_parents)\n",
              "bias_function $(brkga_data.bias)\n",
              "independent_populations $(brkga_data.num_independent_populations)\n",
              # "pr_minimum_distance $(brkga_data.pr_minimum_distance)\n",
              # "pr_type $(brkga_data.pr_type)\n",
              # "alpha_block_size $(brkga_data.alpha_block_size)\n",
              # "pr_percentage $(brkga_data.pr_percentage)\n",
              "exchange_interval $(external_params.exchange_interval)\n",
              "num_exchange_indivuduals $(external_params.num_exchange_indivuduals)\n",
              "reset_interval $(external_params.reset_interval)\n",
        )
    end
    nothing
end
