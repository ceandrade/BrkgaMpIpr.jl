################################################################################
# types_io_tests.jl: unit tests for I/O type handling.
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

@testset "Parsing BiasFunction" begin
    @test parse(BrkgaMpIpr.BiasFunction, "CONSTANT") == CONSTANT
    @test parse(BrkgaMpIpr.BiasFunction, "CUBIC") == CUBIC
    @test parse(BrkgaMpIpr.BiasFunction, "EXPONENTIAL") == EXPONENTIAL
    @test parse(BrkgaMpIpr.BiasFunction, "LINEAR") == LINEAR
    @test parse(BrkgaMpIpr.BiasFunction, "LOGINVERSE") == LOGINVERSE
    @test parse(BrkgaMpIpr.BiasFunction, "QUADRATIC") == QUADRATIC
    @test parse(BrkgaMpIpr.BiasFunction, "QuAdRaTiC") == QUADRATIC
    @test parse(BrkgaMpIpr.BiasFunction, "Custom") == CUSTOM
    @test_throws ArgumentError parse(BrkgaMpIpr.BiasFunction, "invalid")
end

################################################################################

@testset "Parsing PathRelinkingType" begin
    @test parse(BrkgaMpIpr.PathRelinkingType, "DIRECT") == DIRECT
    @test parse(BrkgaMpIpr.PathRelinkingType, "direct") == DIRECT
    @test parse(BrkgaMpIpr.PathRelinkingType, "PERMUTATION") == PERMUTATION
    @test parse(BrkgaMpIpr.PathRelinkingType, "permutation") == PERMUTATION
    @test_throws ArgumentError parse(BrkgaMpIpr.PathRelinkingType, "invalid")
end

################################################################################

@testset "Parsing PathRelinkingSelection" begin
    @test parse(BrkgaMpIpr.PathRelinkingSelection, "BESTSOLUTION") == BESTSOLUTION
    @test parse(BrkgaMpIpr.PathRelinkingSelection, "bestsolution") == BESTSOLUTION
    @test parse(BrkgaMpIpr.PathRelinkingSelection, "RANDOMELITE") == RANDOMELITE
    @test parse(BrkgaMpIpr.PathRelinkingSelection, "randomelite") == RANDOMELITE
    @test_throws ArgumentError parse(BrkgaMpIpr.PathRelinkingSelection, "invalid")
end

################################################################################

@testset "load_configuration()" begin
    config_path = joinpath(@__DIR__, "configuration_files")

    @test_throws LoadError load_configuration(".")
    @test_throws SystemError load_configuration("")
    @test_throws LoadError load_configuration(joinpath(config_path, "missing_value.conf"))
    @test_throws LoadError load_configuration(joinpath(config_path, "unknown_param.conf"))
    @test_throws LoadError load_configuration(joinpath(config_path, "wrong_type.conf"))
    @test_throws LoadError load_configuration(joinpath(config_path, "missing_param.conf"))
    @test_throws LoadError load_configuration(joinpath(config_path, "wrong_bias_function.conf"))
    @test_throws LoadError load_configuration(joinpath(config_path, "wrong_pr_selection.conf"))

    brkga_params, control_params =
        load_configuration(joinpath(config_path, "regular.conf"))

    @test brkga_params.population_size == 500
    @test brkga_params.elite_percentage == 0.30
    @test brkga_params.mutants_percentage == 0.15
    @test brkga_params.num_elite_parents == 2
    @test brkga_params.total_parents == 3
    @test brkga_params.bias_type == LOGINVERSE
    @test brkga_params.num_independent_populations == 3
    @test brkga_params.pr_number_pairs == 0
    @test brkga_params.pr_minimum_distance == 0.15
    @test brkga_params.pr_type == PERMUTATION
    @test brkga_params.pr_selection == RANDOMELITE
    @test brkga_params.alpha_block_size == 1.0
    @test brkga_params.pr_percentage == 1.0
    @test control_params.exchange_interval == 200
    @test control_params.num_exchange_indivuduals == 2
    @test control_params.reset_interval == 600
end

################################################################################

@testset "write_configuration()" begin
    #########################
    # From config file
    #########################

    config_path = joinpath(@__DIR__, "configuration_files")

    brkga_params, external_params = load_configuration(
        joinpath(config_path, "regular.conf")
    )

    @test_throws SystemError write_configuration("/invalid", brkga_params,
                                                 external_params)

    temp_filename = tempname()
    write_configuration(temp_filename, brkga_params, external_params)

    result = ""
    open(temp_filename) do file
        result = lowercase(read(file, String))
    end
    rm(temp_filename)

    standard = """
population_size 500
elite_percentage 0.3
mutants_percentage 0.15
num_elite_parents 2
total_parents 3
bias_type loginverse
num_independent_populations 3
pr_number_pairs 0
pr_minimum_distance 0.15
pr_type permutation
pr_selection randomelite
alpha_block_size 1.0
pr_percentage 1.0
exchange_interval 200
num_exchange_indivuduals 2
reset_interval 600
"""
    @test result == standard

    local_param_values = [
        default_param_values[param_index["instance"]],
        default_param_values[param_index["decode!"]],
        default_param_values[param_index["opt_sense"]],
        default_param_values[param_index["seed"]],
        default_param_values[param_index["chr_size"]]
    ]

    brkga_data, external_params =
        build_brkga(local_param_values...,
                    joinpath(config_path, "regular.conf"))

    #########################
    # From direct building
    #########################
    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)
    brkga_params = brkga_data.params
    external_params = ExternalControlParams(100, 200, 300)

    temp_filename = tempname()
    write_configuration(temp_filename, brkga_params, external_params)

    result = ""
    open(temp_filename) do file
        result = lowercase(read(file, String))
    end
    rm(temp_filename)

    standard = """
population_size 10
elite_percentage 0.3
mutants_percentage 0.1
num_elite_parents 1
total_parents 2
bias_type loginverse
num_independent_populations 3
pr_number_pairs 0
pr_minimum_distance 0.0
pr_type direct
pr_selection bestsolution
alpha_block_size 1.0
pr_percentage 1.0
exchange_interval 100
num_exchange_indivuduals 200
reset_interval 300
"""

    @test result == standard

    #########################
    # From direct building
    #########################
    param_values = deepcopy(default_param_values)
    brkga_data = build_brkga(param_values...)
    set_bias_custom_function!(brkga_data, x -> 1 / x)

    brkga_params = brkga_data.params
    external_params = ExternalControlParams(100, 200, 300)

    temp_filename = tempname()
    write_configuration(temp_filename, brkga_params, external_params)

    result = ""
    open(temp_filename) do file
        result = lowercase(read(file, String))
    end
    rm(temp_filename)

    standard = """
population_size 10
elite_percentage 0.3
mutants_percentage 0.1
num_elite_parents 1
total_parents 2
bias_type custom
num_independent_populations 3
pr_number_pairs 0
pr_minimum_distance 0.0
pr_type direct
pr_selection bestsolution
alpha_block_size 1.0
pr_percentage 1.0
exchange_interval 100
num_exchange_indivuduals 200
reset_interval 300
"""
    @test result == standard
end
