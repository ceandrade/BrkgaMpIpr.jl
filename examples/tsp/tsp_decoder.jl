################################################################################
# tsp_decoder.jl: simple permutation decoder for the Traveling Salesman Problem.
#
# (c) Copyright 2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Dec 28, 2018 by ceandrade
# Last update: Dec 28, 2018 by ceandrade
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

"""
    function tsp_decode!(chromosome::Array{Float64}, instance::Instance,
                         rewrite::Bool = true)::Float64

Simple Traveling Salesman Problem decoder. It creates a permutation of nodes
induced by the chromosome and computes the cost of the tour.

**NOTE** that the method signature is exact the sample required by BRKGA
methods, but the instance type is TSP_Instance.
"""
################################################################################

function tsp_decode!(chromosome::Array{Float64}, instance::TSP_Instance,
                     rewrite::Bool = true)::Float64

    permutation = Array{Tuple{Float64, Int64}}(undef, instance.num_nodes)
    for (index, key) in enumerate(chromosome)
        permutation[index] = (key, index)
    end

    sort!(permutation)

    cost = distance(instance, permutation[1][2], permutation[end][2])
    for i in 1:(instance.num_nodes - 1)
        cost += distance(instance, permutation[i][2], permutation[i + 1][2])
    end

    return cost
end
