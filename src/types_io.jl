################################################################################
# types_io.jl: Input/output/parsing methods for internal data strucutures.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 24, 2018 by ceandrade
# Last update: Dec 27, 2018 by ceandrade
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

Parse `value` returning a valid `BiasFunction` enumeration.

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
    elseif value == "CUSTOM"
        return CUSTOM
    end

    throw(ArgumentError("cannot parse $value as BiasFunction"))
end

################################################################################

"""
    parse(::Type{PathRelinkingType}, value::String)::PathRelinkingType

Parse `value` returning a valid `PathRelinkingType` enumeration.

# Throws
- `ArgumentError`: in case the type description does not match.
"""
function parse(::Type{PathRelinkingType}, value::String)::PathRelinkingType
    value = uppercase(strip(value))
    if value == "DIRECT"
        return DIRECT
    elseif value == "PERMUTATION"
        return PERMUTATION
    end

    throw(ArgumentError("cannot parse $value as PathRelinkingType"))
end

################################################################################

"""
    parse(::Type{PathRelinkingSelection}, value::String)::PathRelinkingSelection

Parse `value` returning a valid `PathRelinkingSelection` enumeration.

# Throws
- `ArgumentError`: in case the selection description does not match.
"""
function parse(::Type{PathRelinkingSelection},
               value::String)::PathRelinkingSelection
    value = uppercase(strip(value))
    if value == "BESTSOLUTION"
        return BESTSOLUTION
    elseif value == "RANDOMELITE"
        return RANDOMELITE
    end

    throw(ArgumentError("cannot parse $value as PathRelinkingSelection"))
end

################################################################################

"""
    write_configuration(filename::String, brkga_data::BrkgaData,
            external_params::ExternalControlParams = ExternalControlParams())

Write the parameters from `brkga_data` and `external_params` into `filename`.

# Throws
- `ArgumentError`: in case the bias description does not match.


"""
function write_configuration(filename::String, brkga_data::BrkgaData,
        external_params::ExternalControlParams = ExternalControlParams())

    elite_percentage = brkga_data.elite_size / brkga_data.population_size;
    mutants_percentage = brkga_data.num_mutants / brkga_data.population_size;

    output_string = """
population_size $(brkga_data.population_size)
elite_percentage $(elite_percentage)
mutants_percentage $(mutants_percentage)
mutants_percentage $(mutants_percentage)
elite_parents $(brkga_data.num_elite_parents)
total_parents $(brkga_data.total_parents)
bias_function $(brkga_data.bias)
independent_populations $(brkga_data.num_independent_populations)
pr_number_pairs $(brkga_data.pr_number_pairs)
pr_minimum_distance $(brkga_data.pr_minimum_distance)
pr_type $(brkga_data.pr_type)
pr_selection $(brkga_data.pr_selection)
alpha_block_size $(brkga_data.alpha_block_size)
pr_percentage $(brkga_data.pr_percentage)
exchange_interval $(external_params.exchange_interval)
num_exchange_indivuduals $(external_params.num_exchange_indivuduals)
reset_interval $(external_params.reset_interval)
"""

    open(filename, "w") do file
        write(file, output_string)
    end
    nothing
end
