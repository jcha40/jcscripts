function maximum_subarray(vec::Vector{<:Number}, min_length::Int64, max_length::Int64)::Tuple{Float64, Int64, Int64}
    """
    Find the maximum sum of a subarray of length between min_length and max_length.
    """
    rolling_sum = 0.
    for i in 1:min_length - 1
        rolling_sum += vec[i]
    end
    curr = rolling_sum
    curr_len = min_length - 1
    max_sum = -Inf
    max_len = 0
    max_idx = 0
    for i in min_length:length(vec)
        next = if curr_len == max_length
            curr + vec[i] - vec[i - max_length]
        else
            curr_len += 1
            curr + vec[i]
        end
        
        rolling_sum += vec[i]
        curr = if next >= rolling_sum
            next
        else
            curr_len = min_length
            rolling_sum
        end
        rolling_sum -= vec[i - min_length + 1]

        if curr > max_sum
            max_sum = curr
            max_len = curr_len
            max_idx = i
        end
    end
    return max_sum, max_idx - max_len + 1, max_idx
end