################################################################################
# types_io.jl: Input/output/parsing methods for internal data strucutures.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 24, 2018 by ceandrade
# Last update: Jan 04, 2018 by ceandrade
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

Parse `value` returning a valid [`BiasFunction`](@ref) enumeration.

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

Parse `value` returning a valid [`PathRelinkingType`](@ref) enumeration.

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

Parse `value` returning a valid [`PathRelinkingSelection`](@ref) enumeration.

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
    load_configuration(configuration_file::String)::
            Tuple{BrkgaParams, ExternalControlParams}

Load the parameters from `configuration_file` returning them as a tuple.

# Throws
- `LoadError`: in cases of the file is an invalid configuration file,
  parameters are missing, or parameters are ill-formatted.
- `SystemError`: in case the configuration files cannot be openned.
"""
function load_configuration(configuration_file::String)::
        Tuple{BrkgaParams, ExternalControlParams}

    # Create a dictionaty with fields and their types.
    param_names_types = Dict(
        [name => typ for (name, typ) in zip(fieldnames(BrkgaParams),
                                            BrkgaParams.types)]
    )
    push!(param_names_types,
        [name => typ for (name, typ)
         in zip(fieldnames(ExternalControlParams),
                ExternalControlParams.types)]...
    )

    param_given = Dict([name => false for name in keys(param_names_types)])

    brkga_params = BrkgaParams()
    control_params = ExternalControlParams()

    lines = Array{String,1}()
    open(configuration_file) do file
        lines = readlines(file)
    end
    if length(lines) == 0
        throw(LoadError(configuration_file, 0,
                        "cannot read '$configuration_file'"))
    end

    for (line_number, line) in enumerate(lines)
        line = strip(line)
        if length(line) == 0 || line[1] == '#'
            continue
        end

        param_name = ""
        value = 0
        try
            param_name, value = split(line)
            param_name = lowercase(param_name)

            field = Symbol(param_name)
            if field in fieldnames(BrkgaParams)
                data = brkga_params
            elseif field in fieldnames(ExternalControlParams)
                data = control_params
            else
                throw(KeyError(""))
            end

            setfield!(data, field, parse(param_names_types[field],
                                         String(value)))
            param_given[field] = true
        catch err
            if isa(err, BoundsError)
                throw(LoadError(configuration_file, line_number, "error line " *
                                "$line_number of '$configuration_file': " *
                                "missing parameter or value"))

            elseif isa(err, KeyError)
                throw(LoadError(configuration_file, line_number,
                                "parameter '$param_name' unknown"))

            elseif isa(err, ArgumentError)
                throw(LoadError(configuration_file, line_number,
                                "invalid value for '$param_name': $value"))
            else
                throw(err)
            end
        end
    end

    missing_params = ""
    for (name, value) in param_given
        if !value
            missing_params *= "'$name',"
        end
    end
    if length(missing_params) > 0
        throw(LoadError(configuration_file, 0,
                        "missing parameters: $missing_params"))
    end

    return (brkga_params, control_params)
end

################################################################################

"""
    function write_configuration(filename::String, brkga_params::BrkgaParams,
                                 external_params::ExternalControlParams)

Write `brkga_params` and `external_params` into `filename`.

# Throws
- `SystemError`: in case the configuration files cannot be openned.
"""
function write_configuration(filename::String, brkga_params::BrkgaParams,
                             external_params::ExternalControlParams)

    output_string = ""
    for field in fieldnames(BrkgaParams)
        output_string *= "$field $(getfield(brkga_params, field))\n"
    end
    for field in fieldnames(ExternalControlParams)
        output_string *= "$field $(getfield(external_params, field))\n"
    end

    open(filename, "w") do file
        write(file, output_string)
    end
    nothing
end

"""
    function write_configuration(filename::String, brkga_data::BrkgaData,
                                 external_params::ExternalControlParams =
                                                        ExternalControlParams())

Write the parameters from `brkga_data.params` and `external_params`
into `filename`.

# Throws
- `SystemError`: in case the configuration files cannot be openned.
"""
function write_configuration(filename::String, brkga_data::BrkgaData,
                             external_params::ExternalControlParams =
                                                        ExternalControlParams())
    write_configuration(filename, brkga_data.params, external_params)
    nothing
end
