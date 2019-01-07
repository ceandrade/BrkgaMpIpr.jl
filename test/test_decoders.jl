################################################################################
# TestDecoders.jl: decoders for tests.
#
# (c) Copyright 2019, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on:  Apr 20, 2018 by ceandrade
# Last update: Dec 18, 2018 by ceandrade
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

function sum_decode!(chromosome::Array{Float64}, instance::Instance,
                 rewrite::Bool = true)::Float64
    tmp = chromosome + instance.data
    tmp /= maximum(tmp)
    if rewrite
        chromosome .= tmp
    end
    return Float64(sum(tmp))
end

################################################################################

function rank_decode!(chromosome::Array{Float64}, instance::Instance,
                      rewrite::Bool = true)::Float64
    tmp = chromosome + instance.data

    rank = 0
    @inbounds for i in 1:(length(tmp) - 1)
        if tmp[i] < tmp[i + 1]
            rank += 1
        end
    end

    if rewrite
        chromosome .= tmp
    end
    return Float64(rank)
end