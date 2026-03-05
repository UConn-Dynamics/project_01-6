module Simulate

using OrdinaryDiffEq
using ..Physics

function pendulum!(du, u, p, t)
    θ  = u[1]
    θd = u[2]

    L, w1, h1, g, Ω = p
    
    du[1] = θd
    du[2] = θdd_numeric(θ, θd, L, w1, h1, g, Ω, t)
end

function run_sim(IC, params, tspan)
    prob = ODEProblem(pendulum!, IC, tspan, params)
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8)
    return sol
end

export run_sim

end