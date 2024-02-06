using SpecialFunctions

function logbetah(a::Float64, b::Float64; n::Int64 = 20, loghalf::Float64 = log(.5))::Float64
    """
    Compute the logarithm of the regularized half-beta function using the continued fraction representation.
    """
    if a < b
        return log1p(-exp(logbetah(b, a, n=n)))
    end
    frac = 1.
    for k in n:-1:1
        frac = 1 - (((a + k - 1) * (a + b + k - 1) * .5) / ((a + k * 2 - 2) * (a + k * 2 - 1))) / (1 + ((k * (b - k) * .5) / ((a + k * 2 - 1) * (a + k * 2))) / frac)
    end
    return a * loghalf + b * loghalf - log(a) - logbeta(a, b) - log(frac)
end