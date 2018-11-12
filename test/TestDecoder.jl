function decode!(chromosome::Array{Float64}, instance::Instance,
                 rewrite::Bool = true)::Float64
    tmp = chromosome + instance.data
    tmp /= maximum(tmp)
    if rewrite
        # **NOTE:** "chromosome = copy(tmp)" is too slow.
        # Tt is better do it manually.
        @inbounds for i in 1:length(tmp)
            chromosome[i] = tmp[i]
        end
    end
    return sum(tmp)
end
