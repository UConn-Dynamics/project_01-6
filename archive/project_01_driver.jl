using OrdinaryDiffEq
using Symbolics
using Plots
# ------------------------------------------------------------------------
# define physical constants 
# ------------------------------------------------------------------------
const g  = 9.81     
const h1 = 0.20
const w1 = 0.10
const L  = 0.15

# ------------------------------------------------------------------------
# define rotation speeds
# ------------------------------------------------------------------------
const Ω_slow = 1.0     # rad/s
const Ω_fast = 12.0    # rad/s

# ------------------------------------------------------------------------
# define the ODE (first-order system)
# ------------------------------------------------------------------------
function pendulum!(du, u, p, t)
    θ  = u[1]
    θd = u[2]

    Ω = p

    du[1] = θd
    du[2] = (Ω^2 / L) * (w1 + L * sin(θ)) * cos(θ) - (g / L) * sin(θ)
end

# ------------------------------------------------------------------------
# UPDATE CODE TO DO HAND DERIVATION 
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# define wrapper function for pendulum simulation 
# ------------------------------------------------------------------------

function run_sim(θ_0, θd_0, Ω, tspan(0.0, 10.0))
    u_0 = [θ_0, θd_0]
    prob = ODEProblem(pendulum!, u_0, tspan, Ω)
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8)
    return sol
end

# ------------------------------------------------------------------------
# record positions from a solution
# ------------------------------------------------------------------------

function xyz_from_sol(sol, Ω)
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


# ------------------------------------------------------------------------
# post processing
# ------------------------------------------------------------------------
function post_process(sol, Ω)
    t, x, y, z = xyz_from_sol(sol, Ω)
    plt1 = plot()
    plot(t, x, label="x(t)", xlabel="Time (s)", ylabel="Position (m)")
    plot!(t, y, label="y(t)")
    plot!(t, z, label="z(t)")
end


# ------------------------------------------------------------------------
# define initial conditions
# ------------------------------------------------------------------------
function define_initial_conditions()
    θ_0 = 0.0
    θd_0 = 0.0
    return θ_0, θd_0
end

# ------------------------------------------------------------------------
# main driver (for now)
# ------------------------------------------------------------------------
function main()
    #
end

@main()
