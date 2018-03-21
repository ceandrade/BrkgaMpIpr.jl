function evolve!(brkga_data::BrkgaData, problem_instance::AbstractInstance,
                 decode!::Function)

    # results = Array{Float64}(brkga_data.num_chromosomes)
    #
    # Threads.@threads for i in eachindex(brkga_data.chromosomes)
    #     results[i] = decode!(brkga_data.chromosomes[i], problem_instance)
    # end
end
