module TestDecoder
export decode!

import TestInstance

function decode!(chromosome::Array{Float64}, instance::TestInstance.Instance,
                 rewrite::Bool = false)::Float64

    chromosome[:] = sqrt.(chromosome + instance.data)

    total = 0
    last = chromosome[1]
    for i = 1:size(chromosome)[1]
        if last < chromosome[i]
            total += 1
        end
        last = chromosome[i]
    end
    return Float64(total)
end

end
