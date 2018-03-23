################################################################################
# BrkgaMpIpr.jl: Main module for BRKGA-MP-IPR framework.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Mar 23, 2018 by ceandrade
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

module BrkgaMpIpr

include("types.jl")
include("init.jl")
include("evolution.jl")

# TODO (ceandrade): exporting of each item of the enums is weird.
# If Julia change these in the future, remove the items from exporting.
export Sense, MINIMIZE, MAXIMIZE

# TODO (ceandrade): exporting of each item of the enums is weird.
# If Julia change these in the future, remove the items from exporting.
export BiasFunction, CONSTANT, CUBIC, EXPONENTIAL, LINEAR, LOGINVERSE, QUADRATIC

export BrkgaData, AbstractInstance, ExternalControlParams
export build_brkga, set_bias_custom_function!
export evolve!

end
