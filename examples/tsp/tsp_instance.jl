################################################################################
# tsp_instance.jl: data structures and support function to deal with instances
# of the Traveling Salesman Problem.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
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
    struct TSP_Instance

Represents an instance for the Traveling Salesman Problem. The constructor
loads a upper triangular matrix as following:

    number of nodes (n)
    dist12 dist13 dist14 ... dist1n
    dist23 dist24 ... dist2(n - 1)
    ...
    dist(n-2)(n-1)

For example, for n = 4 we have

    4
    12 13 14
    23 24
    34

**NOTE** that this structure inherits from `AbstractInstance` conforming
required signatures.
"""
struct TSP_Instance <: AbstractInstance
    num_nodes::Int64
    distances::Array{Float64}

    function TSP_Instance(filename::String)
        lines = Array{String,1}()
        open(filename) do file
            lines = readlines(file)
        end
        if length(lines) == 0
            throw(LoadError(filename, 0, "cannot read '$filename'"))
        end

        num_nodes = 0
        distances = Array{Float64, 1}()
        line_number = 0
        try
            num_nodes = parse(Int64, lines[1])

            matrix_size::Int64 = (num_nodes * (num_nodes - 1)) / 2
            distances = Array{Float64, 1}(undef, matrix_size)
            for i in 1:(num_nodes - 1)
                line_number = i + 1
                values = split(lines[i + 1])
                base::Int64 = ((i - 1) * num_nodes) - ((i - 1) * i / 2)
                for j in (i + 1):num_nodes
                    distances[base + (j - i)] = parse(Float64, values[j - i])
                end
            end
        catch err
            throw(LoadError(filename, line_number,
                            "error line $(line_number) of '$filename'"))
        end

        new(num_nodes, distances)
    end
end

################################################################################

"""
    function distance(instance::TSP_Instance, i::Int64, j::Int64)

Return the distance between nodes `i` and `j`.
"""
@inline function distance(instance::TSP_Instance, i::Int64, j::Int64)
    if i > j
        i, j = j, i
    end
    return instance.distances[((i - 1) * instance.num_nodes) -
                              Int64(((i - 1) * i / 2)) + (j - i)]
end
