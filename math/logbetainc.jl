using SpecialFunctions

function logbetainc(a::Real, b::Real, x::Real; n::Int64 = 20, min_x_swap::Float64 = 1e-16)::Float64
    """
    Compute the logarithm of the regularized incomplete beta function using the continued fraction representation.
    """
    if x < 0. || x > 1. || a <= 0. || b <= 0.
        return NaN
    end
    # Using the identity I_x(a, b) = 1 - I_{1-x}(b, a) converges faster when x > (a + 1) / (a + b + 2)
    # However, this can cause underflow when x is very small, so we only swap when x >= min_x_swap
    if ((a + 1.) * (1. - x) < (b + 1.) * x) && x >= min_x_swap
        return log1p(-exp(logbetainc(b, a, 1. - x, n=n)))
    end
    frac = 1.
    for k in n:-1:1
        frac = 1 - (((a + k - 1) * (a + b + k - 1) * x) / ((a + k * 2 - 2) * (a + k * 2 - 1))) / (1 + ((k * (b - k) * x) / ((a + k * 2 - 1) * (a + k * 2))) / frac)
    end
    return a * log(x) + b * log1p(-x) - log(a) - logbeta(a, b) - log(frac)
end

function logbetah(a::Real, b::Real; n::Int64 = 20, loghalf::Float64 = log(.5))::Float64
    """
    Compute the logarithm of the regularized half-beta function using the continued fraction representation.
    """
    if a <= 0. || b <= 0.
        return NaN
    end
    if a < b
        return log1p(-exp(logbetah(b, a, n=n)))
    end
    frac = 1.
    for k in n:-1:1
        frac = 1 - (((a + k - 1) * (a + b + k - 1) * .5) / ((a + k * 2 - 2) * (a + k * 2 - 1))) / (1 + ((k * (b - k) * .5) / ((a + k * 2 - 1) * (a + k * 2))) / frac)
    end
    return a * loghalf + b * loghalf - log(a) - logbeta(a, b) - log(frac)
end