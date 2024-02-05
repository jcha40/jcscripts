using SpecialFunctions

function logbetai(x::Float64, a::Float64, b::Float64; n::Int64 = 20)::Float64
    """
    Compute the logarithm of the incomplete regularized beta function using the continued fraction representation.
    """
    if a < b
        return log1p(-exp(logbetai(x, b, a, n=n)))
    end
    frac = 1.
    for k in n:-1:1
        frac = 1 - (((a + k - 1) * (a + b + k - 1) * x) / ((a + k * 2 - 2) * (a + k * 2 - 1))) / (1 + ((k * (b - k) * x) / ((a + k * 2 - 1) * (a + k * 2))) / frac)
    end
    return a * log(x) + b * log(1 - x) - log(a) - logbeta(a, b) - log(frac)
end