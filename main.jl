include("src/derive_equations.jl")
include("src/physics.jl")
include("src/simulate.jl")
include("src/visualize.jl")

using .DeriveEquations
using .Physics
using .Simulate
using .Visualize
using Plots

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

function main()
    # physical constants
    L = 0.15   # length of pendulum
    w1 = 0.10  # horizontal offset of pivot
    h1 = 0.20  # vertical offset of pivot
    g = 9.81   # gravitational acceleration
    tspan = (0.0, 10.0)  # time span for simulation

    Ω_slow = 1.0     # rad/s
    Ω_fast = 12.0    # rad/s

    # initial conditions
    IC_zero = [0.0, 0.0]    # initial angle of zero
    IC_small = [0.05, 0.0]  # small initial angle
    
    params_slow = (L, w1, h1, g, Ω_slow)
    params_fast = (L, w1, h1, g, Ω_fast)

    # run simulations
    sol_slow_zero = run_sim(IC_zero, params_slow, tspan)
    sol_slow_small = run_sim(IC_small, params_slow, tspan)
    sol_fast_zero = run_sim(IC_zero, params_fast, tspan)
    sol_fast_small = run_sim(IC_small, params_fast, tspan)

    # case automation
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

    # plot results
    #θ(t) plots
    for (cases, filename, title_label) in [(cases_theta_slow, "results/theta_vs_time_low_omega.png", "Low Ω"),
                                           (cases_theta_fast, "results/theta_vs_time_high_omega.png", "High Ω")]
        p_theta = plot(title="θ(t) - $title_label", xlabel="t (s)", ylabel="θ (rad)", legend=:outertopright)
        for (label, sol) in cases
            plot!(p_theta, sol.t, getindex.(sol.u,1), label=label, lw=2)
        end
        display(p_theta)
        savefig(p_theta, filename)
        println("Press Enter to continue to next plot...")
        readline()
    end

    # 3D trajectories
    for (cases, filename, title_label) in [(cases_3d_slow, "results/3d_trajectories_low_omega.png", "Low Ω"),
                                           (cases_3d_fast, "results/3d_trajectories_high_omega.png", "High Ω")]
        p_3d = plot3d(title="3D Trajectories - $title_label", xlabel="x", ylabel="y", zlabel="z", legend=:outertopright)
        for (label, sol, Ω) in cases
            t, x, y, z = xyz_from_sol(sol, L, w1, h1, Ω)
            plot3d!(x, y, z, label=label, lw=2)
        end
        display(p_3d)
        savefig(p_3d, filename)
        println("Press Enter to continue to next plot...")
        readline()
    end

    # Δθ''(t) plots
    p_delta = plot(title="Δθ''(t) Between Symbolic Function and Hand-Calculated EOM", xlabel="t (s)", ylabel="Δθ'' (rad/s²)")
    for (label, sol, Ω) in cases_delta
        ts, deltas = solution_delta(sol, L, w1, h1, g, Ω)
        plot!(p_delta, ts, deltas, label=label, lw=2)
    end
    display(p_delta)
    savefig(p_delta, "results/theta_double_dot_delta.png")
    println("Press Enter to exit...")
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

