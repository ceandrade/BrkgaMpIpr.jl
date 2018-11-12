import Random

mutable struct Instance <: AbstractInstance
    data::Array{Float64}

    function Instance(size::Int64)
        rng = Random.MersenneTwister(1234)
        new(rand(rng, size))
    end
end
