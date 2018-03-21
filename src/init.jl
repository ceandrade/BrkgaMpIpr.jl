################################################################################
# init.jl: Initialization routines.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Mar 20, 2018 by ceandrade
# Last update: Mar 21, 2018 by ceandrade
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

function init(opt_sense::Sense, num_chromosomes::Int64,
              chromosome_size::Int64)::BrkgaData

    brkga_data = BrkgaData()
    brkga_data.opt_sense = opt_sense

    brkga_data.rng = MersenneTwister(1234)

    # brkga_data.num_chromosomes = num_chromosomes
    # brkga_data.chromosome_size = chromosome_size
    #
    # brkga_data.chromosomes = Array{Array{Float64}}(0)
    # for i in 1:num_chromosomes
    #     push!(brkga_data.chromosomes, rand(brkga_data.rng, chromosome_size))
    # end

    return brkga_data
end
