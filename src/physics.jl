module Physics

using ModelingToolkit
using ..DeriveEquations: t, θ, θdd_rhs, θdd_func

# ------------------------------------------------------------------------
# Symbolic Parameters
# ------------------------------------------------------------------------

@parameters L w1 h1 g Ω

θt = θ

# ------------------------------------------------------------------------
# Symbolic Expression of Hand-Calculated θ''
# ------------------------------------------------------------------------

manual_expr = (Ω^2 / L) * (w1 + L * sin(θt)) * cos(θt) - (g / L) * sin(θt)

# ------------------------------------------------------------------------
# Numeric Functions
# ------------------------------------------------------------------------

# numeric θ'' from symbolic function
θdd_numeric(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val) = 
    θdd_func(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val)

# hand calculation for θ''
θdd_manual_numeric(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val) = 
    (Ω_val^2 / L_val) * (w1_val + L_val * sin(θ_val)) * cos(θ_val) - (g_val / L_val) * sin(θ_val)

# difference between symbolic function and hand calculation
θdd_delta_numeric(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val) =
    θdd_numeric(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val) - 
    θdd_manual_numeric(θ_val, θd_val, L_val, w1_val, h1_val, g_val, Ω_val, t_val)

# ------------------------------------------------------------------------
# Exported Functions
# ------------------------------------------------------------------------

export θdd_numeric, θdd_manual_numeric, θdd_delta_numeric

end
