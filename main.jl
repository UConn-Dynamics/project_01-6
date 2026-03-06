# ------------------------------------------------------------------------
# Main Driver Script for Pendulum Simulation
# ------------------------------------------------------------------------
# Includes modules, runs simulations, generates plots/animations, 
# computes deltas
# ------------------------------------------------------------------------

include("src/derive_equations.jl")
include("src/physics.jl")
include("src/simulate.jl")
include("src/visualize.jl")

using .DeriveEquations
using .Physics
using .Simulate
using .Visualize
using Plots

# ------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------

"""
    solution_delta(sol, L, w1, h1, g, Ω)

Computes the delta between the symbolic and hand-calculated EOMs

# Arguments
- `sol`: solution object containing 't' and 'u' arrays (from `run_sim`)
- `L, w1, h1, g, Ω`: pendulum parameters

# Returns
- `(ts, deltas)`: time array and corresponding delta values

"""

function solution_delta(sol, L, w1, h1, g, Ω)
    ts = sol.t
    θs = getindex.(sol.u, 1)
    θds = getindex.(sol.u, 2)

    deltas = [
        θdd_delta_numeric(θs[i], θds[i], L, w1, h1, g, Ω, ts[i])
        for i in eachindex(ts)
    ]

    return ts, deltas
end

"""
    safe_filename(label::String)

Converts a label string into a safe filename (converts '.' to '')    
"""

function safe_filename(label::String)
    omega = match(r"Ω=([\d.]+)", label).captures[1] |> x -> replace(x, "." => "")
    theta = match(r"θ₀=([\d.]+)", label).captures[1] |> x -> replace(x, "." => "")
    return "omega_$(omega)_theta0_$(theta)"
end

# ------------------------------------------------------------------------
# Main Simulation and Visualization
# ------------------------------------------------------------------------

function main()
    # --------------------------------------------------------------------
    # physical constants
    # --------------------------------------------------------------------
    
    L = 0.15   # length of pendulum
    w1 = 0.10  # horizontal offset of pivot
    h1 = 0.20  # vertical offset of pivot
    g = 9.81   # gravitational acceleration
    tspan = (0.0, 10.0)  # time span for simulation

    Ω_slow = 1.0     # rad/s
    Ω_fast = 12.0    # rad/s

    # --------------------------------------------------------------------
    # Initial Conditions
    # --------------------------------------------------------------------

    IC_zero = [0.0, 0.0]    # initial angle of zero
    IC_small = [0.25, 0.0]  # small initial angle
    
    params_slow = (L, w1, h1, g, Ω_slow)
    params_fast = (L, w1, h1, g, Ω_fast)

    # --------------------------------------------------------------------
    # Run Simulations
    # --------------------------------------------------------------------

    sol_slow_zero = run_sim(IC_zero, params_slow, tspan)
    sol_slow_small = run_sim(IC_small, params_slow, tspan)
    sol_fast_zero = run_sim(IC_zero, params_fast, tspan)
    sol_fast_small = run_sim(IC_small, params_fast, tspan)

    # --------------------------------------------------------------------
    # Organize Cases for Plots and Animations
    # --------------------------------------------------------------------

    cases_theta_slow = [
        ("Ω=1.0, θ₀=0.0", sol_slow_zero),
        ("Ω=1.0, θ₀=0.05", sol_slow_small)
    ]

    cases_theta_fast = [
        ("Ω=12.0, θ₀=0.0", sol_fast_zero),
        ("Ω=12.0, θ₀=0.05", sol_fast_small)
    ]

    cases_3d_slow = [
        ("Ω=1.0, θ₀=0.0", sol_slow_zero, Ω_slow),
        ("Ω=1.0, θ₀=0.05", sol_slow_small, Ω_slow)
    ]

    cases_3d_fast = [
        ("Ω=12.0, θ₀=0.0", sol_fast_zero, Ω_fast),
        ("Ω=12.0, θ₀=0.05", sol_fast_small, Ω_fast)
    ]

    cases_delta = [
        ("Ω=1.0, θ₀=0.0", sol_slow_zero, Ω_slow),
        ("Ω=1.0, θ₀=0.05", sol_slow_small, Ω_slow),
        ("Ω=12.0, θ₀=0.0", sol_fast_zero, Ω_fast),
        ("Ω=12.0, θ₀=0.05", sol_fast_small, Ω_fast)
    ]

    # --------------------------------------------------------------------
    # Plotting Limits
    # --------------------------------------------------------------------

    xlims = (-0.27, 0.27)
    ylims = (-0.27, 0.27)
    zlims = (0.0, 0.37)

    # --------------------------------------------------------------------
    # Static Plots
    # --------------------------------------------------------------------

    #θ(t) plots
    for (cases, filename, title_label) in [(cases_theta_slow, "results/theta_vs_time_low_omega.png", "Low Ω"),
                                           (cases_theta_fast, "results/theta_vs_time_high_omega.png", "High Ω")]
        p_theta = plot(title="θ(t) - $title_label", xlabel="t (s)", ylabel="θ (rad)", legend=:outertopright)
        for (label, sol) in cases
            plot!(p_theta, sol.t, getindex.(sol.u,1), label=label, lw=2)
        end
        display(p_theta)
        savefig(p_theta, filename)
        #println("Press Enter to continue to next plot...")
        #readline()
    end

    # 3D trajectories
    for (cases, filename, title_label) in [(cases_3d_slow, "results/3d_trajectories_low_omega.png", "Low Ω"),
                                           (cases_3d_fast, "results/3d_trajectories_high_omega.png", "High Ω")]
        p_3d = plot3d(title="3D Trajectories - $title_label", xlabel="x", ylabel="y", zlabel="z", legend=:outertopright, xlims=xlims, ylims=ylims, zlims=zlims)
        for (label, sol, Ω) in cases
            t, x, y, z = xyz_from_sol(sol, L, w1, h1, Ω)
            plot3d!(x, y, z, label=label, lw=2)
        end
        display(p_3d)
        savefig(p_3d, filename)
        #println("Press Enter to continue to next plot...")
        #readline()
    end

    # Δθ''(t) plots
    p_delta = plot(title="Δθ''(t) Between Symbolic Function and Hand-Calculated EOM", xlabel="t (s)", ylabel="Δθ'' (rad/s²)")
    for (label, sol, Ω) in cases_delta
        ts, deltas = solution_delta(sol, L, w1, h1, g, Ω)
        plot!(p_delta, ts, deltas, label=label, lw=2)
    end
    display(p_delta)
    savefig(p_delta, "results/theta_double_dot_delta.png")
    #println("Press Enter to contine to animations...")
    #readline()

    # --------------------------------------------------------------------
    # Animations
    # --------------------------------------------------------------------

    for (label, sol, Ω) in cases_delta
        fname_base = safe_filename(label)

        # animate θ(t)
        animate_theta(sol, Ω, filename="results/theta_$(fname_base).gif")
        #println("Press Enter to contine to next animation...")
        #readline()

        # animate θ'(t)
        animate_thetadot(sol, Ω, filename="results/thetadot_$(fname_base).gif")
        #println("Press Enter to contine to next animation...")
        #readline()


        # animate 3D trajectory
        animate_3d(sol, L, w1, h1, Ω, filename="results/3d_trajectory_$(fname_base).gif")
        #println("Press Enter to contine to next animation...")
        #readline()

        # animate pendulum dashboard
        animate_pendulum_dashboard(sol, L, w1, h1, Ω, filename="results/dashboard_$(fname_base).gif")
        #println("Press Enter to exit...")
        #readline()

    end

end

# ------------------------------------------------------------------------
# Main Driver Call
# ------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

