using SpecialFunctions

function logbetai(x::Float64, a::Float64, b::Float64; n::Int64 = 100)::Float64
    frac, s = if n % 2 == 0
        1, n
    else
        k = div(n, 2)
        1 - ((a + k) * (a + b + k) * x) / ((a + n - 1) * (a + n)), n - 1
    end
    for i in s:-2:2
        k = div(i, 2)
        frac = 1 - (((a + k - 1) * (a + b + k - 1) * x) / ((a + i - 2) * (a + i - 1))) / (1 + ((k * (b - k) * x) / ((a + i - 1) * (a + i))) / frac)
    end
    return a * log(x) + b * log(1 - x) - log(a) - logbeta(a, b) - log(frac)
end