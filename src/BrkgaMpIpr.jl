module BrkgaMpIpr
export AbstractInstance, BrkgaData, init, evolve!

abstract type AbstractInstance end

mutable struct BrkgaData
    num_chromosomes::Int64
    chromosome_size::Int64
    chromosomes::Array{Array{Float64}}
    rng::MersenneTwister
    BrkgaData() = new()
end

function init(num_chromosomes::Int64, chromosome_size::Int64)
    brkga_data = BrkgaData()
    brkga_data.rng = MersenneTwister(1234)

    brkga_data.num_chromosomes = num_chromosomes
    brkga_data.chromosome_size = chromosome_size

    brkga_data.chromosomes = Array{Array{Float64}}(0)
    for i in 1:num_chromosomes
        push!(brkga_data.chromosomes, rand(brkga_data.rng, chromosome_size))
    end

    return brkga_data
end

function evolve!(brkga_data::BrkgaData, problem_instance::AbstractInstance,
                 decode!::Function)

    results = Array{Float64}(brkga_data.num_chromosomes)

    Threads.@threads for i in eachindex(brkga_data.chromosomes)
        results[i] = decode!(brkga_data.chromosomes[i], problem_instance)
    end
end

end
