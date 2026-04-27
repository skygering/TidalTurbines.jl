module TidalTurbines

# Write your package code here.
# -------------------------------------------------------------------
# Turbine definitions
# -------------------------------------------------------------------
const turbines = [
    (
        x0 = 3.80,
        y0 = 1.70,
        D = 0.20,
        r_support = 0.15,
        Ct_rated = 0.80,
        Cp_rated = 0.45,
        u_in = 0.20,
        u_rated = 0.80,
        u_out = 2.50,
        h_min = 0.02
    )
    # ,
    # (
    #     x0 = 4.20,
    #     y0 = 1.70,
    #     D = 0.20,
    #     r_support = 0.15,
    #     Ct_rated = 0.80,
    #     Cp_rated = 0.45,
    #     u_in = 0.20,
    #     u_rated = 0.80,
    #     u_out = 2.50,
    #     h_min = 0.02
    # )
]

const XI_BUMP = 1.45661

# -------------------------------------------------------------------
# Smooth step function
# -------------------------------------------------------------------
@inline function smoothstep(s)
    if s <= 0
        return 0.0
    elseif s >= 1
        return 1.0
    else
        return s^2 * (3 - 2s)
    end
end

# -------------------------------------------------------------------
# Bump function for turbine footprint
# -------------------------------------------------------------------
@inline function bump_1d(x, p, r)
    ξ = (x - p) / r
    if abs(ξ) < 1
        return exp(1 - 1 / (1 - ξ^2))
    else
        return zero(x)
    end
end

@inline function turbine_density_single(x, turb)
    x1, x2 = x
    ψx = bump_1d(x1, turb.x0, turb.r_support)
    ψy = bump_1d(x2, turb.y0, turb.r_support)
    return ψx * ψy / (XI_BUMP * turb.r_support^2)
end

# -------------------------------------------------------------------
# Turbine thrust and power curves
# -------------------------------------------------------------------
@inline function Ct_turbine(U, turb)
    u_in     = turb.u_in
    u_rated  = turb.u_rated
    u_out    = turb.u_out
    Ct_rated = turb.Ct_rated

    if U <= u_in
        return 0.0
    elseif U < u_rated
        s = (U - u_in) / (u_rated - u_in)
        return Ct_rated * smoothstep(s)
    elseif U < u_out
        s = (U - u_rated) / (u_out - u_rated)
        return Ct_rated * (1.0 - 0.5 * smoothstep(s))
    else
        return 0.0
    end
end

@inline function Cp_turbine(U, turb)
    u_in     = turb.u_in
    u_rated  = turb.u_rated
    u_out    = turb.u_out
    Cp_rated = turb.Cp_rated

    if U <= u_in
        return 0.0
    elseif U < u_rated
        s = (U - u_in) / (u_rated - u_in)
        return Cp_rated * smoothstep(s)
    elseif U < u_out
        return Cp_rated * (u_rated / U)^3
    else
        return 0.0
    end
end

# -------------------------------------------------------------------
# Combined source terms:
#   Manning friction + distributed turbine drag
# -------------------------------------------------------------------
@inline function source_terms_manning_plus_turbines(u, x, t,
                                                    equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    T = eltype(x)
    zero_T = zero(T)

    # -------------------------------
    # Manning friction
    # -------------------------------
    n = 0.001
    h_eff = max(h, equations.threshold_limiter)

    h_fric = (h_eff^2 + max(h_eff^2, 1e-8)) / (2 * h_eff)
    Sf = -equations.gravity * n^2 * h_fric^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    S_h   = zero_T
    S_hv1 = Sf * hv_1
    S_hv2 = Sf * hv_2
    S_b   = zero_T

    # -------------------------------
    # Turbine drag from all turbines
    # -------------------------------
    if h_eff > equations.threshold_limiter
        v1 = hv_1 / h_eff
        v2 = hv_2 / h_eff
        U  = sqrt(v1^2 + v2^2)

        for turb in turbines
            if h_eff > turb.h_min
                dA = turbine_density_single(x, turb)
                if dA > 0
                    At = 0.25 * pi * turb.D^2
                    Ct = Ct_turbine(U, turb)

                    Kt = 0.5 * Ct * At * dA / h_eff

                    S_hv1 -= Kt * U * hv_1
                    S_hv2 -= Kt * U * hv_2
                end
            end
        end
    end

    return SVector(S_h, S_hv1, S_hv2, S_b)
end

# -------------------------------------------------------------------
# Turbine power density diagnostics
# -------------------------------------------------------------------
@inline function turbine_power_density_single(u, x, equations, turb;
                                              rho = 1000.0)
    h, hv_1, hv_2, _ = u

    h_eff = max(h, equations.threshold_limiter)
    if h_eff <= turb.h_min
        return 0.0
    end

    dA = turbine_density_single(x, turb)
    if dA <= 0
        return 0.0
    end

    v1 = hv_1 / h_eff
    v2 = hv_2 / h_eff
    U  = sqrt(v1^2 + v2^2)

    Cp = Cp_turbine(U, turb)
    At = 0.25 * pi * turb.D^2

    return 0.5 * rho * Cp * At * dA * U^3
end

@inline function turbine_power_density_total(u, x, equations; rho = 1000.0)
    p = 0.0
    for turb in turbines
        p += turbine_power_density_single(u, x, equations, turb; rho = rho)
    end
    return p
end

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_terms_manning_plus_turbines);

                                    using DelimitedFiles
using StaticArrays
using OrdinaryDiffEq

# ============================================================
# Simple turbine power output using nearest DG node to center
# ============================================================

# ---- extract one nodal conservative state
@inline function get_node_state(u_wrap, i, j, element)
    return SVector(u_wrap[1, i, j, element],
                   u_wrap[2, i, j, element],
                   u_wrap[3, i, j, element],
                   u_wrap[4, i, j, element])
end

# ---- get node physical coordinates
# NOTE:
# This works for the common Trixi cache layout used by DGSEM.
# If your local Trixi version uses a slightly different field name,
# send me `fieldnames(typeof(semi.cache.elements))` and I can adjust it.
@inline function get_node_xy(semi, i, j, element)
    x = semi.cache.elements.node_coordinates[1, i, j, element]
    y = semi.cache.elements.node_coordinates[2, i, j, element]
    return SVector(x, y)
end

# ---- find nearest DG node for each turbine center
function build_turbine_probe_map(semi, turbines)
    solver = semi.solver
    mesh   = semi.mesh

    probe_map = Vector{NamedTuple}(undef, length(turbines))

    for n in eachindex(turbines)
        turb = turbines[n]
        best_dist2 = Inf
        best_entry = (element = 0, i = 0, j = 0)

        for element in eachelement(mesh)
            for j in eachnode(solver), i in eachnode(solver)
                x = get_node_xy(semi, i, j, element)
                dx = x[1] - turb.x0
                dy = x[2] - turb.y0
                dist2 = dx^2 + dy^2
                if dist2 < best_dist2
                    best_dist2 = dist2
                    best_entry = (element = element, i = i, j = j)
                end
            end
        end

        probe_map[n] = best_entry
    end

    return probe_map
end

# ---- compute per-turbine powers from nearest-node velocities
function compute_turbine_powers_simple(semi, u_ode, probe_map; rho = 1000.0)
    u_wrap = Trixi.wrap_array(u_ode, semi)
    eqs    = semi.equations

    powers = zeros(Float64, length(turbines))

    for n in eachindex(turbines)
        turb = turbines[n]
        probe = probe_map[n]

        u_node = get_node_state(u_wrap, probe.i, probe.j, probe.element)
        h, hv_1, hv_2, _ = u_node

        h_eff = max(h, eqs.threshold_limiter)

        if h_eff <= turb.h_min
            powers[n] = 0.0
            continue
        end

        v1 = hv_1 / h_eff
        v2 = hv_2 / h_eff
        U  = sqrt(v1^2 + v2^2)

        Cp = Cp_turbine(U, turb)
        At = 0.25 * pi * turb.D^2

        # simple actuator-disk style turbine power estimate
        powers[n] = 0.5 * rho * Cp * At * U^3
    end

    return powers, sum(powers)
end

# ---- logger object
mutable struct SimpleTurbinePowerLogger
    semi
    probe_map
    filename::String
    dt::Float64
    next_time::Float64
    last_time::Float64
    last_power::Float64
    cumulative_energy::Float64
    rho::Float64
    initialized::Bool
end

function SimpleTurbinePowerLogger(semi, turbines;
                                  filename = "out/turbine_power.csv",
                                  dt = 0.1,
                                  rho = 1000.0)
    probe_map = build_turbine_probe_map(semi, turbines)
    return SimpleTurbinePowerLogger(semi, probe_map, filename, dt,
                                    0.0, 0.0, 0.0, 0.0, rho, false)
end

# ---- initialize csv
function initialize_simple_power_logger!(logger::SimpleTurbinePowerLogger)
    mkpath(dirname(logger.filename))

    header = ["time"]
    append!(header, ["P_turbine_$k" for k in 1:length(turbines)])
    push!(header, "P_total")
    push!(header, "E_total")

    open(logger.filename, "w") do io
        writedlm(io, permutedims(header), ',')
    end

    logger.initialized = true
end

# ---- write one row
function write_simple_power_row!(logger::SimpleTurbinePowerLogger, t, powers, P_total)
    row = vcat([t], powers, [P_total, logger.cumulative_energy])
    open(logger.filename, "a") do io
        writedlm(io, permutedims(row), ',')
    end
end

# ---- callback action
function simple_power_logger_affect!(integrator, logger::SimpleTurbinePowerLogger)
    if !logger.initialized
        initialize_simple_power_logger!(logger)
    end

    t = integrator.t
    powers, P_total = compute_turbine_powers_simple(logger.semi, integrator.u,
                                                    logger.probe_map; rho = logger.rho)

    # trapezoidal cumulative energy
    if t > logger.last_time
        logger.cumulative_energy += 0.5 * (logger.last_power + P_total) * (t - logger.last_time)
    end

    write_simple_power_row!(logger, t, powers, P_total)

    logger.last_time = t
    logger.last_power = P_total
    logger.next_time += logger.dt
end

# ---- callback constructor
function SimplePowerOutputCallback(semi, turbines;
                                   filename = "out/turbine_power.csv",
                                   dt = 0.1,
                                   rho = 1000.0)

    logger = SimpleTurbinePowerLogger(semi, turbines;
                                      filename = filename,
                                      dt = dt,
                                      rho = rho)

    condition = (u, t, integrator) -> t >= logger.next_time
    affect!   = integrator -> simple_power_logger_affect!(integrator, logger)

    return DiscreteCallback(condition, affect!)
end

# create callback here
power_callback = SimplePowerOutputCallback(semi, turbines;
                                           filename = "out/turbine_power.csv",
                                           dt = 0.05,
                                           rho = 1000.0)

callbacks = CallbackSet(analysis_callback,
                        stepsize_callback,
                        save_solution,
                        power_callback);
using Trixi
using TrixiShallowWater
using TrixiBottomTopography

include("bathymetry.jl")
include("boundary_conditions.jl")
include("source_terms.jl")

end
