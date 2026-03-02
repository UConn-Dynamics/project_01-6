module Visualize

using Plots

function xyz_from_sol(sol, L, w1, h1, Ω)
    t = sol.t
    θs = getindex.(sol.u, 1)

    x = similar(t)
    y = similar(t)
    z = similar(t)

    for (i , (ti, θ)) in enumerate(zip(t, θs))
        R = w1 + L * sin(θ)
        x[i] = R * cos(Ω * ti)
        y[i] = R * sin(Ω * ti)
        z[i] = h1 - L * cos(θ)
    end

    return t, x, y, z
end

function plot_θ(sol; label="")
    plot(sol.t, getindex.(sol.u,1), xlabel="t", ylabel="θ", label=label)
end

export xyz_from_sol, plot_θ

end