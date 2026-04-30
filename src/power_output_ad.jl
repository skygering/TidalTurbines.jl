function compute_turbine_powers(semi, u_ode, turbines; rho = 1000.0)
    solver    = semi.solver
    equations = semi.equations
    cache     = semi.cache

    u_wrap = Trixi.wrap_array(u_ode, semi)

    nturb = length(turbines)
    powers = zeros(Float64, nturb)

    for element in eachelement(semi.solver, semi.cache)
        for j in eachnode(solver), i in eachnode(solver)
            x = get_node_xy(cache, i, j, element)
            dV = get_node_volume_weight(solver, cache, i, j, element)
            u_node = get_node_state(u_wrap, i, j, element)

            for n in 1:nturb
                powers[n] += turbine_power_density_single(u_node, x, equations, turbines[n];
                                                          rho = rho) * dV
            end
        end
    end

    return powers, sum(powers)
end

function compute_total_turbine_energy(sol, semi, turbines; rho = 1000.0)
    nt = length(sol.t)
    P_hist = zeros(Float64, nt)

    for k in 1:nt
        _, P_total = compute_turbine_powers(semi, sol.u[k], turbines; rho = rho)
        P_hist[k] = P_total
    end

    E = 0.0
    for k in 2:nt
        dt = sol.t[k] - sol.t[k - 1]
        E += 0.5 * (P_hist[k] + P_hist[k - 1]) * dt
    end

    return E, P_hist
end

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

# @inline function get_node_xy(semi, i, j, element)
#     x = semi.cache.elements.node_coordinates[1, i, j, element]
#     y = semi.cache.elements.node_coordinates[2, i, j, element]
#     return SVector(x, y)
# end

@inline function get_node_xy(cache, i, j, element)
    x = cache.elements.node_coordinates[1, i, j, element]
    y = cache.elements.node_coordinates[2, i, j, element]
    return SVector(x, y)
end

@inline function get_node_volume_weight(solver, cache, i, j, element)
    w_i = solver.basis.weights[i]
    w_j = solver.basis.weights[j]
    J   = cache.elements.inverse_jacobian[i, j, element] # may need adaptation
    return w_i * w_j / J
end

# ---- find nearest DG node for each turbine center
# function build_turbine_probe_map(semi, turbines)
#     solver = semi.solver
#     mesh   = semi.mesh

#     probe_map = Vector{NamedTuple}(undef, length(turbines))

#     for n in eachindex(turbines)
#         turb = turbines[n]
#         best_dist2 = Inf
#         best_entry = (element = 0, i = 0, j = 0)

#         for element in eachelement(mesh)
#             for j in eachnode(solver), i in eachnode(solver)
#                 x = get_node_xy(semi, i, j, element)
#                 dx = x[1] - turb.x0
#                 dy = x[2] - turb.y0
#                 dist2 = dx^2 + dy^2
#                 if dist2 < best_dist2
#                     best_dist2 = dist2
#                     best_entry = (element = element, i = i, j = j)
#                 end
#             end
#         end

#         probe_map[n] = best_entry
#     end

#     return probe_map
# end

# function compute_turbine_powers(semi, u_ode, turbines; rho = 1000.0)
#     solver    = semi.solver
#     equations = semi.equations
#     cache     = semi.cache
#     mesh      = semi.mesh

#     # Wrap ODE vector into nodal array used by Trixi
#     u_wrap = Trixi.wrap_array(u_ode, semi)

#     nturb = length(turbines)
#     powers = zeros(Float64, nturb)

#     # for element in eachelement(mesh)
#     for element in eachelement(semi.solver, semi.cache)
#         # println("Computing turbine power for element")
#         for j in eachnode(solver), i in eachnode(solver)
#             # println("  node (i=$i, j=$j) ...",i,j)
#             x = get_node_xy(cache, i, j, element)
#             dV = get_node_volume_weight(solver, cache, i, j, element)
#             u_node = get_node_state(u_wrap, i, j, element)

#             for n in 1:nturb
#                 powers[n] += turbine_power_density_single(u_node, x, equations, turbines[n];
#                                                           rho = rho) * dV
#                 # println("    turbine $n power contribution: ", powers[n])
#             end
#         end
#     end

#     return powers, sum(powers)
# end

# ---- compute per-turbine powers from nearest-node velocities
# function compute_turbine_powers_simple(semi, u_ode, probe_map, turbines; rho=1000.0)
#     u_wrap = Trixi.wrap_array(u_ode, semi)
#     eqs    = semi.equations

#     powers = zeros(Float64, length(turbines))

#     for n in eachindex(turbines)
#         turb = turbines[n]
#         probe = probe_map[n]

#         u_node = get_node_state(u_wrap, probe.i, probe.j, probe.element)
#         h, hv_1, hv_2, _ = u_node

#         h_eff = max(h, eqs.threshold_limiter)

#         if h_eff <= turb.h_min
#             powers[n] = 0.0
#             continue
#         end

#         v1 = hv_1 / h_eff
#         v2 = hv_2 / h_eff
#         U  = sqrt(v1^2 + v2^2)

#         Cp = Cp_turbine(U, turb)
#         At = 0.25 * pi * turb.D^2

#         # simple actuator-disk style turbine power estimate
#         powers[n] = 0.5 * rho * Cp * At * U^3
#     end

#     return powers, sum(powers)
# end

# ---- logger object
mutable struct SimpleTurbinePowerLogger
    semi
    turbines
    # probe_map
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
    # probe_map = build_turbine_probe_map(semi, turbines)
    return SimpleTurbinePowerLogger(semi, turbines, filename, dt,
                                    0.0, 0.0, 0.0, 0.0, rho, false)
end

# ---- initialize csv
function initialize_simple_power_logger!(logger::SimpleTurbinePowerLogger)
    mkpath(dirname(logger.filename))

    header = ["time"]
    append!(header, ["P_turbine_$k" for k in 1:length(logger.turbines)])
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
    powers, P_total = compute_turbine_powers(logger.semi, integrator.u,
                                                    logger.turbines;
                                                    rho = logger.rho)

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
# power_callback = SimplePowerOutputCallback(semi, turbines;
#                                            filename = "out/turbine_power.csv",
#                                            dt = 0.05,
#                                            rho = 1000.0)

# callbacks = CallbackSet(analysis_callback,
#                         stepsize_callback,
#                         save_solution,
#                         power_callback);