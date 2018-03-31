################################################################################
# types_io_tests.jl: unit tests for I/O type handling.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 24, 2018 by ceandrade
# Last update: Mar 31, 2018 by ceandrade
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
    @test_throws ArgumentError parse(BrkgaMpIpr.BiasFunction, "invalid")
end

################################################################################

@testset "write_configuration()" begin
    #########################
    # From config file
    #########################

    config_path = joinpath(@__DIR__, "configuration_files")

    const local_param_values = [
        default_param_values[param_index["instance"]],
        default_param_values[param_index["decode!"]],
        default_param_values[param_index["opt_sense"]],
        default_param_values[param_index["seed"]],
        default_param_values[param_index["chr_size"]]
    ]

    brkga_data, external_params =
        build_brkga(local_param_values...,
                    joinpath(config_path, "regular.conf"))

    @test_throws SystemError write_configuration(brkga_data, external_params,
                                                 "/invalid")

    temp_filename = tempname()
    write_configuration(brkga_data, external_params, temp_filename)

    result = ""
    open(temp_filename) do file
        result = lowercase(readstring(file))
    end
    rm(temp_filename)

    # TODO (ceandrade): add the path relink parameters.
    standard = """
population_size 500
elite_percentage 0.3
mutants_percentage 0.15
mutants_percentage 0.15
elite_parents 2
total_parents 3
bias_function loginverse
independent_populations 3
exchange_interval 200
num_exchange_indivuduals 2
reset_interval 600
"""
    @test result == standard

    #########################
    # From direct building
    #########################
    param_values = copy(default_param_values)
    brkga_data = build_brkga(param_values...)

    temp_filename = tempname()
    write_configuration(brkga_data, ExternalControlParams(0, 0, 0),
                        temp_filename)

    result = ""
    open(temp_filename) do file
        result = lowercase(readstring(file))
    end
    rm(temp_filename)

    standard = """
population_size 10
elite_percentage 0.3
mutants_percentage 0.1
mutants_percentage 0.1
elite_parents 1
total_parents 2
bias_function loginverse
independent_populations 3
exchange_interval 0
num_exchange_indivuduals 0
reset_interval 0
"""
    @test result == standard
end
