using Plots

function bezier(control_points::Tuple{Float64, Float64}...)::Function
    n = length(control_points) - 1
    coefs = ones(Int64, n + 1)
    for i in 1:n
        coefs[i + 1] = div(coefs[i] * (n - i + 1), i)
    end
    function pbezier(t::Float64)::Tuple{Float64, Float64}
        x = 0.
        y = 0.
        for i in 1:n + 1
            b = coefs[i] * ((1 - t) ^ (n - i + 1)) * (t ^ (i - 1))
            x += b * control_points[i][1]
            y += b * control_points[i][2]
        end
        return x, y
    end
    return pbezier
end

paths = [
    "M 0 100 L 33 0 L 66 0 L 100 100 L 75 100 L 66 75 L 33 75 L 25 100 L 0 100 M 41 55 L 58 55 L 50 25 L 41 55",
    "M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 10 75 28 Z",
    "M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 100 48 L 55 48 L 55 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 5 75 28 Z",
    "M 0 0 L 0 20 L 35 20 L 35 100 L 65 100 L 65 20 L 100 20 L 100 0 Z"
]
colors = ["#FF0000", "#0000FF", "#FFA500", "#228B22"]

pattern = r"(?:(?:M|L) (-?\d+ -?\d+))|(?:C (-?\d+ -?\d+ -?\d+ -?\d+ -?\d+ -?\d+))"
dna_symbols = Shape[]
for path in paths
    xy = Tuple{Float64, Float64}[]
    prev = nothing
    for m in eachmatch(pattern, path)
        if m[1] !== nothing
            x, y = parse.(Float64, split(m[1]))
            prev = (x, 100. - y)
            push!(xy, prev)
        else
            x1, y1, x2, y2, x3, y3 = parse.(Float64, split(m[2]))
            f = bezier(prev, (x1, 100. - y1), (x2, 100. - y2), (x3, 100. - y3))
            push!(xy, f.((1:100) ./ 100.)...)
            prev = (x3, 100. - y3)
        end
    end
    push!(dna_symbols, Shape(xy))
end

function plot_logo(pwm::Matrix{Float64}, alpha::Float64 = 0.)
    p = plot(legend=false, axis=false, grid=false, aspect_ratio=:equal)
    pwm_width = size(pwm)[1]
    for x in 1:pwm_width
        pwv = pwm[x, :]
        n_sites = round(sum(pwv))
        vec = (pwv .+ alpha) ./ (n_sites + 4. * alpha)
        vec[vec .> 0] .*= log2.(vec[vec .> 0])
        bits = 2. + sum(vec)
        heights = pwv .* (bits / 2. * 200. / n_sites)
        idx = sortperm(pwv)

        y_offset = 0.
        for i in idx
            plot!(p, translate(Plots.scale(dna_symbols[i], 1., heights[i] / 100., (50., 0.)), x * 100., y_offset), fillcolor=colors[i], linecolor=nothing)
            y_offset += heights[i]
        end
    end
    return p
end

function plot_logo_stack(aligned_pwms::Array{Float64, 3})
    p = plot(legend=false, axis=false, grid=false, aspect_ratio=:equal)
    n_pwms, pwm_width, _ = size(aligned_pwms)
    y = 0.
    for i in 1:n_pwms
        pwm = aligned_pwms[i, :, :]
        for x in 1:pwm_width
            pwv = pwm[x, :]
            if !isnan(pwv[1])
                n_sites = round(sum(pwv))
                vec = pwv ./ n_sites
                vec[vec .> 0] .*= log2.(vec[vec .> 0])
                bits = 2. + sum(vec)
                heights = pwv .* (bits / 2. * 200. / n_sites)
                idx = sortperm(pwv)

                y_offset = 0.
                for j in idx
                    plot!(p, translate(Plots.scale(dna_symbols[j], 1., heights[j] / 100., (50., 0.)), x * 100., y + y_offset), fillcolor=colors[j], linecolor=nothing)
                    y_offset += heights[j]
                end
            end
        end
        y -= 210.
    end
    return p
end