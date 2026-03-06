module Visualize

using Plots
using Measures

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

function build_frame(t_i, mass_x_i, mass_y_i, mass_z_i, w1, h1, Ω)
    pivot_x = w1 * cos(Ω * t_i)
    pivot_y = w1 * sin(Ω * t_i)
    pivot_z = h1

    # rod
    rod_x = [pivot_x, mass_x_i]
    rod_y = [pivot_y, mass_y_i]
    rod_z = [pivot_z, mass_z_i] 

    #frame
    frame_x = [0.0, 0.0]
    frame_y = [0.0, 0.0]
    frame_z = [0.0, h1]

    # horizontal support
    arm_x = [0.0, pivot_x]
    arm_y = [0.0, pivot_y]
    arm_z = [h1, h1]

    return rod_x, rod_y, rod_z, pivot_x, pivot_y, pivot_z, frame_x, frame_y, frame_z, arm_x, arm_y, arm_z

end

function plot_θ(sol; label="")
    plot(sol.t, getindex.(sol.u,1), xlabel="t", ylabel="θ", label=label)

end
# ------------------------------------
# Individual plot animations
# ------------------------------------

# animate θ(t)
function animate_theta(sol, Ω; filename="results/theta.gif", fps=30)
    θs = getindex.(sol.u, 1)
    θmin  = minimum(θs)
    θmax  = maximum(θs)

    #inds = 1:10:length(sol.t)
    inds = 1:length(sol.t)
    anim = @animate for i in inds
        t = sol.t[1:i]
        θs_i = θs[1:i]
        p = plot(t, θs_i, xlabel="t (s)", ylabel="θ (rad)", title="θ(t): Ω=$(Ω), θ₀ =$(θs_i[1])", legend=false, ylim=(θmin-0.05, θmax+0.05))
        scatter!(p, [t[end]], [θs_i[end]], color=:red, left_margin=5mm, right_margin=5mm, bottom_margin=5mm, top_margin=10mm, size = (600, 400))
        
    end

    gif(anim, filename, fps=fps)
    #display(anim)

end

# animate θ'(t)
function animate_thetadot(sol, Ω; filename="results/thetadot.gif", fps=30)
    θs  = getindex.(sol.u, 1)
    θd_s = getindex.(sol.u, 2)
    θdmin  = minimum(θd_s)
    θdmax  = maximum(θd_s)
   
    #inds = 1:10:length(sol.t)
    inds = 1:length(sol.t)
    anim = @animate for i in inds
        t = sol.t[1:i]
        θds = θd_s[1:i]
        p = plot(t, θds, xlabel="t (s)", ylabel="θ'", title="θ'(t): Ω=$(Ω), θ₀ =$(θs[1])", legend=false, ylim=(θdmin-0.05, θdmax+0.05))
        scatter!(p, [t[end]], [θds[end]], color=:blue, left_margin=5mm, right_margin=5mm, bottom_margin=5mm, top_margin=10mm, size = (600, 400))
        
    end

    gif(anim, filename, fps=fps)
    #display(anim)

end

function animate_3d(sol, L, w1, h1, Ω; filename="results/3d_trajectory.gif", fps=30)
    t, x, y, z = xyz_from_sol(sol, L, w1, h1, Ω)
    θs = getindex.(sol.u, 1)

    #inds = 1:10:length(sol.t)
    inds = 1:length(sol.t)
    anim = @animate for i in inds
        
        t_i = t[i]
        mass_x_i, mass_y_i, mass_z_i = x[i], y[i], z[i]
        
        # build frame for visualization
        rod_x, rod_y, rod_z, pivot_x, pivot_y, pivot_z, frame_x, frame_y, frame_z, arm_x, arm_y, arm_z = build_frame(t_i, mass_x_i, mass_y_i, mass_z_i, w1, h1, Ω)

        p = plot3d(x[1:i], y[1:i], z[1:i], xlabel="x", ylabel="y", zlabel="z", title="3D Trajectory: Ω=$(Ω), θ₀ =$(θs[1])", label="Position Trace", lw = 2, color=:green, linestyle = :dot, left_margin=5mm, right_margin=5mm, bottom_margin=5mm, top_margin=10mm)
        scatter3d!(p, [mass_x_i], [mass_y_i], [mass_z_i], color=:red, label="Mass", markersize=5)
        scatter3d!(p,  [pivot_x], [pivot_y], [pivot_z], color=:black, label="Pivot", markersize=5)

        plot3d!(p, rod_x, rod_y, rod_z, color=:blue, label="Rod", lw=3)
        plot3d!(p, frame_x, frame_y, frame_z, color=:black, label="Vertical Support", lw=3)
        plot3d!(p, arm_x, arm_y, arm_z, color=:gray, label="Horizontal Support", lw=3)
        plot3d!(p, legend=:outertopright)
        #plot3d!(p, title="3D Trajectory", legend=:outertopright, xlims=(-L-w1-0.5, L+w1+0.5), ylims=(-L-w1-0.5, L+w1+0.5), zlims=(0, h1+L+0.5))
        
    end

    gif(anim, filename, fps=fps)
    #display(anim)

end

# combined dashboard
function animate_pendulum_dashboard(sol, L, w1, h1, Ω; filename="results/pendulum_dashboard.gif", fps=30)
    θs = getindex.(sol.u, 1)
    θmin  = minimum(θs)
    θmax  = maximum(θs)

    θd_s = getindex.(sol.u, 2)
    θdmin  = minimum(θd_s)
    θdmax  = maximum(θd_s)

    t, x, y, z = xyz_from_sol(sol, L, w1, h1, Ω)

    #inds = 1:10:length(sol.t)
    inds = 1:length(sol.t)
    anim = @animate for i in inds
        
        θs_i  = θs[1:i]
        θds_i = θd_s[1:i]
        t_i = t[i]
        mass_x_i, mass_y_i, mass_z_i = x[i], y[i], z[i]

        # θ plot
        p1 = plot(t[1:i], θs_i, xlabel="t (s)", ylabel="θ (rad)", title="θ(t)", legend=false, ylim=(θmin-0.05, θmax+0.05), left_margin=10mm, right_margin=5mm, bottom_margin=5mm, top_margin=5mm)
        scatter!(p1, [t_i], [θs_i[end]], color=:red)

        # θ' plot
        p2 = plot(t[1:i], θds_i, xlabel="t (s)", ylabel="θ' (rad/s)", title="θ'(t)", legend=false, ylim=(θdmin-0.05, θdmax+0.05), left_margin=10mm, right_margin=5mm, bottom_margin=5mm, top_margin=5mm)
        scatter!(p2, [t_i], [θds_i[end]], color=:blue)

        # build frame for visualization
        rod_x, rod_y, rod_z, pivot_x, pivot_y, pivot_z, frame_x, frame_y, frame_z, arm_x, arm_y, arm_z = build_frame(t_i, mass_x_i, mass_y_i, mass_z_i, w1, h1, Ω)

        p3 = plot(x[1:i], y[1:i], z[1:i], lw=3, color=:green, linestyle = :dot, title="3D Trajectory", label="Position Trace", left_margin=10mm, right_margin=5mm, bottom_margin=5mm, top_margin=5mm)
        scatter3d!(p3, [mass_x_i], [mass_y_i], [mass_z_i], color=:red, label="Mass", markersize=5)
        scatter3d!(p3,  [pivot_x], [pivot_y], [pivot_z], color=:black, label="Pivot", markersize=5)

        plot3d!(p3, rod_x, rod_y, rod_z, color=:blue, label="Rod", lw=3)
        plot3d!(p3, frame_x, frame_y, frame_z, color=:black, label="Vertical Support", lw=3)
        plot3d!(p3, arm_x, arm_y, arm_z, color=:black, label="Horizontal Support", lw=3)
        
        plot3d!(p3, xlabel="x", ylabel="y", zlabel="z", legend=:outertopright)
        #plot3d!(p3, xlabel="x", ylabel="y", zlabel="z", xlims=(-L-w1-0.5, L+w1+0.5), ylims=(-L-w1-0.5, L+w1+0.5), zlims=(0, h1+L+0.5), legend=:outertopright)
        plot(p1, p2, p3, layout=(3,1), size=(800,1200), suptitle="Pendulum Simulation Dashboard: Ω=$(Ω), θ₀ =$(θs[1])", top_margin=5mm)
        
    end

    gif(anim, filename, fps=fps)
    #display(anim)

end

export xyz_from_sol, build_frame, plot_θ, animate_theta, animate_thetadot, animate_3d, animate_pendulum_dashboard

end