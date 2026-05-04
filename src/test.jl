
@inline function Cp_from_Ct(Ct)
    Ct_clamped = clamp(Ct, 0.0, 1.0)
    return 0.5 * (1.0 + sqrt(1.0 - Ct_clamped)) * Ct_clamped
end

const Ct_rated_val = 0.516
const Cp_rated_val = Cp_from_Ct(Ct_rated_val)

const turbines = [
    (
        x0 = 3.80,
        y0 = 1.70,
        D = 0.3,
        r_support = 0.5,
        Ct_rated = Ct_rated_val,
        Cp_rated = Cp_rated_val,
        u_in = 1.00,
        u_rated = 3.05,
        u_out = 6.00,
        h_min = 0.02
    )
]

@inline function smoothstep(s)
    if s <= 0
        return 0.0
    elseif s >= 1
        return 1.0
    else
        return s^2 * (3 - 2s)
    end
end

@inline function Ct_turbine(U, turb)
    u_in      = turb.u_in       # 1.0
    u_rated   = turb.u_rated    # 3.05
    u_out     = turb.u_out
    Ct_rated  = turb.Ct_rated   # 0.516

    if U <= 0.0
        return 0.0

    elseif U < u_in
        return Ct_rated * exp(5.0 * (U / u_in - 1.0))

    elseif U <= u_rated
        return Ct_rated

    elseif U < u_out
        Cp = Cp_capped_above_rated(U, turb)
        return Ct_from_Cp_cardano(Cp)
    else
        return 0.0
    end
end

@inline function Cp_turbine(U, turb)
    Ct = Ct_turbine(U, turb)
    return Cp_from_Ct(Ct)
end

@inline function Cp_capped_above_rated(U, turb)
    return turb.Cp_rated * (turb.u_rated / U)^3
end

@inline function Ct_from_Cp_cardano(Cp)
    Cp_clamped = clamp(Cp, 0.0, 16.0 / 27.0)
    θ = acos(1.0 - 27.0 * Cp_clamped / 8.0) / 3.0
    a = (2.0 / 3.0) * (1.0 - cos(θ))
    return 4.0 * a * (1.0 - a)
end

U = 0:0.01:5.5
turb = turbines[1]
Ct = Ct_turbine.(U, Ref(turb))

using Plots
pl1 = Plots.plot(U, Ct, xlabel="U", ylabel="Ct", label="Ct(U)", lineWidth=2)
Cp = Cp_turbine.(U, Ref(turb))

@inline function power_turbine(U, turb)
    A = π * turb.D^2 / 4
    Cp = Cp_turbine(U, turb)
    return 0.5 * 1025.0 * A * Cp * U^3 / 1e6 # MW
end

P = [power_turbine(u, turb) for u in U]
pl2 = Plots.plot(U, Cp, xlabel="U", ylabel="Power [MW]", label="P(U)", 
    legend=:topleft, lineWidth=2)

display(pl1)
display(pl2)