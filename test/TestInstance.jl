module TestInstance
export Instance

import BrkgaMpIpr.AbstractInstance

mutable struct Instance <: AbstractInstance
    data::Array{Float64}

    function Instance(size::Int64)
        rng = MersenneTwister(1234)
        new(rand(rng, size))
    end
end
end
