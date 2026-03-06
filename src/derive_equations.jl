module DeriveEquations

using ModelingToolkit

# ------------------------------------------------------------------------
# Symbolic Variables and Parameters
# ------------------------------------------------------------------------

@variables t θ(t)
@parameters m L w1 h1 g Ω

Dt = Differential(t)
Dθ = Differential(θ)
Dθd = Differential(Dt(θ))

# ------------------------------------------------------------------------
# Coordinate Transformation
# ------------------------------------------------------------------------

# rotation matrix about z
R(ϕ) = [ cos(ϕ) -sin(ϕ) 0;
         sin(ϕ)  cos(ϕ) 0;
         0       0      1]

# local coordinates of mass in rotating frame
r_local = [w1 + L*sin(θ);
           0;
           h1 - L*cos(θ)]

# global coordinates of mass in inertial frame
r = R(Ω*t) * r_local

#velocity of mass in inertial frame
v = expand_derivatives.(Dt.(r))

# ------------------------------------------------------------------------
# Energies and Lagrangian
# ------------------------------------------------------------------------

# kinetic energy
T = simplify((1/2) * m * sum(v .* v))

# potential energy
V = m * g * (h1 - L*cos(θ))

# Lagrangian
Lagr = T - V

# ------------------------------------------------------------------------
# Equations of Motion
# ------------------------------------------------------------------------

# Euler-Lagrange equation
EL_expr = simplify(
    expand_derivatives(
        Dt(Dθd(Lagr)) - Dθ(Lagr)
    )
)

EL_eq = EL_expr ~ 0

# solve for θ''(t)
θdd_rhs = solve_for(EL_eq, Dt(Dt(θ)))[1] |> simplify

# build numeric function for ODE solver
θdd_func = build_function(
    θdd_rhs,
    θ, Dt(θ), L, w1, h1, g, Ω, t;
    expression=Val{false}
)

# ------------------------------------------------------------------------
# Exported Functions
# ------------------------------------------------------------------------

export θdd_func, θdd_rhs, Lagr

end 