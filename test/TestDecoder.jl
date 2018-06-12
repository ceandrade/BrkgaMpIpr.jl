module TestDecoder
export decode!

import TestInstance

function decode!(chromosome::Array{Float64}, instance::TestInstance.Instance,
                 rewrite::Bool = true)::Float64

    tmp = chromosome + instance.data
    tmp /= maximum(tmp)

    if rewrite
        # **NOTE:** the following comment assignments are too slow; it is better
        # do it manually: chromosome = copy(tmp)
        @inbounds for i in 1:length(tmp)
            chromosome[i] = tmp[i]
        end
    end

    return sum(tmp)
end

end
