################################################################################
# types_io_tests.jl: unit tests for I/O type handling.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 24, 2018 by ceandrade
# Last update: Mar 24, 2018 by ceandrade
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
using Base.Test

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
