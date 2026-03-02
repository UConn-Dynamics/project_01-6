include("src/derive_equations.jl")
include("src/physics.jl")
include("src/simulate.jl")
include("src/visualize.jl")

using .DeriveEquations
using .Physics
using .Simulate
using .Visualize
using Plots

function main()
    # physical constants
    L = 0.15   # length of pendulum
    w1 = 0.10  # horizontal offset of pivot
    h1 = 0.20  # vertical offset of pivot
    g = 9.81   # gravitational acceleration

    Ω_slow = 1.0     # rad/s
    Ω_fast = 12.0    # rad/s

    # initial conditions
    IC_zero = [0.0, 0.0]    # initial angle of zero
    IC_small = [0.05, 0.0]  # small initial angle
    
    params_slow = (L, w1, h1, g, Ω_slow)
    params_fast = (L, w1, h1, g, Ω_fast)

    # run simulations
    sol_slow_zero = run_sim(IC_zero, params_slow)
    sol_slow_small = run_sim(IC_small, params_slow)
    sol_fast_zero = run_sim(IC_zero, params_fast)
    sol_fast_small = run_sim(IC_small, params_fast)

    # plot results
    p = plot_θ(sol_slow_zero, label="Ω=1.0, θ₀=0.0")
    plot!(p, sol_slow_small.t, getindex.(sol_slow_small.u,1), label="Ω=1.0, θ₀=0.05")
    plot!(p, sol_fast_zero.t, getindex.(sol_fast_zero.u,1), label="Ω=12.0, θ₀=0.0")
    plot!(p, sol_fast_small.t, getindex.(sol_fast_small.u,1), label="Ω=12.0, θ₀=0.05")
    display(p)
    println("Press Enter to continue to 3D trajectory plots...")
    readline()

    # 3D trajectories
    t_slow_zero, x_slow_zero, y_slow_zero, z_slow_zero = xyz_from_sol(sol_slow_zero, L, w1, h1, Ω_slow)
    t_slow_small, x_slow_small, y_slow_small, z_slow_small = xyz_from_sol(sol_slow_small, L, w1, h1, Ω_slow)
    t_fast_zero, x_fast_zero, y_fast_zero, z_fast_zero = xyz_from_sol(sol_fast_zero, L, w1, h1, Ω_fast)
    t_fast_small, x_fast_small, y_fast_small, z_fast_small = xyz_from_sol(sol_fast_small, L, w1, h1, Ω_fast)
    p3 = plot3d(x_slow_zero, y_slow_zero, z_slow_zero, label="Ω=1.0, θ₀=0.0")
    plot3d!(x_slow_small, y_slow_small, z_slow_small, label="Ω=1.0, θ₀=0.05")
    plot3d!(x_fast_zero, y_fast_zero, z_fast_zero, label="Ω=12.0, θ₀=0.0")
    plot3d!(x_fast_small, y_fast_small, z_fast_small, label="Ω=12.0, θ₀=0.05")
    display(p3)
    println("Press Enter to exit...")
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

