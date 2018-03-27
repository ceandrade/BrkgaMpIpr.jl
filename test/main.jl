################################################################################
# test_main.jl: anchor point for test and development.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Mar 27, 2018 by ceandrade
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

push!(LOAD_PATH, ".")

import TestInstance
import TestDecoder
using BrkgaMpIpr

chr_size = 100
pop_size = 10
elite_percentage = 0.3
mutants_percentage = 0.1
evolutionary_mechanism_on = true
num_elite_parents = 1
total_parents = 2
bias = BrkgaMpIpr.LOGINVERSE
num_independent_populations = 2

print("\n\n>> Building instance")
instance = TestInstance.Instance(chr_size)

print("\n\n>> Building BRKGA data")
# brkga_data = BrkgaMpIpr.build_brkga(instance, TestDecoder.decode!,
#                              BrkgaMpIpr.MAXIMIZE, 2700001, chr_size, pop_size,
#                              elite_percentage, mutants_percentage,
#                              evolutionary_mechanism_on, num_elite_parents,
#                              total_parents, bias, num_independent_populations)

(brkga_data, control_params) =
    BrkgaMpIpr.build_brkga(instance, TestDecoder.decode!,
                           BrkgaMpIpr.MAXIMIZE, 2700001, chr_size,
                           "configuration_files/regular.conf")

print("\n\n>> Initializing BRKGA")
BrkgaMpIpr.initialize!(brkga_data)

print("\n\n>> Reseting")
BrkgaMpIpr.reset!(brkga_data)

# print("\n\n>> Evolving")
# @time BrkgaMpIpr.evolve!(brkga_data, instance, TestDecoder.decode!)

print("\n\n>> Exchanging individuals")
BrkgaMpIpr.exchange_elite!(brkga_data, 2)
