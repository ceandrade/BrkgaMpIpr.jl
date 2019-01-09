################################################################################
# greedy_tour.jl: Simple greedy TSP tour linking the closest nodes in sequence.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Jan 08, 2018 by ceandrade
# Last update: Jan 08, 2018 by ceandrade
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
    function greedy_tour(instance::TSP_Instance)::Array{Int64, 1}

Build a greedy Traveling Salesman Problem tour starting from node 1.

# Returns
The tour cost and a permutation of the nodes representing it.
"""
function greedy_tour(instance::TSP_Instance)::Tuple{Float64, Array{Int64, 1}}
    tour = zeros(Int64, instance.num_nodes)
    tour[1] = 1
    remaining_nodes = Set(collect(2:instance.num_nodes))

    cost = 0.0
    current = 1
    next = 0
    idx = 2
    while length(remaining_nodes) > 0
        best_dist = Inf
        @inbounds for j in remaining_nodes
            dist = distance(instance, current, j)
            if dist < best_dist
                best_dist = dist
                next = j
            end
        end

        cost += best_dist
        tour[idx] = next
        delete!(remaining_nodes, next)
        current = next
        idx += 1
    end
    cost += distance(instance, tour[1], tour[end])

    return (cost, tour)
end
