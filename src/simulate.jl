module Simulate

using OrdinaryDiffEq
using ..Physics

# ------------------------------------------------------------------------
# Pendulum Dynamics
# ------------------------------------------------------------------------

"""
    pendulum!(du, u, p, t)

Defines the ODE system for the pendulum dynamics

# Arguments
- `du`: derivative vector
- `u`: state vector
- `p`: parameters `[L, w1, h1, g, Ω]`
    - `L`: pendulum length
    - `w1, h1`: pivot offsets
    - `g`: gravitational acceleration
    - `Ω`: angular velocity of the rotating frame
- `t`: time

# Modifies
- `du[1]` : θ'
- `du[2]` : θ''

"""

function pendulum!(du, u, p, t)
    θ  = u[1]
    θd = u[2]

    L, w1, h1, g, Ω = p
    
    du[1] = θd
    du[2] = θdd_numeric(θ, θd, L, w1, h1, g, Ω, t)
end

# ------------------------------------------------------------------------
# Simulation Function
# ------------------------------------------------------------------------

"""
    run_sim(IC, params, tspan)

Run the pendulum simulation using `OrdinaryDiffEq`

# Arguments
- `IC`: initial conditions `[θ₀, θ'₀]`
- `params`: parameters `[L, w1, h1, g, Ω]`
- `tspan`: time span `(t_start, t_end)`

# Returns
- `sol`: solution object

"""

function run_sim(IC, params, tspan)
    prob = ODEProblem(pendulum!, IC, tspan, params)
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8)
    return sol
end

# ------------------------------------------------------------------------
# Exported Functions
# ------------------------------------------------------------------------

export run_sim

end