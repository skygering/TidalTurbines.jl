# export initial_condition_tidal_surge
# # ---------------------------
# # Tidal forcing parameters
# # ---------------------------
# const T_tide = 12.42 * 3600
# const ω      = 2π / T_tide
# const A_left   = 1.5 

# # ---------------------------
# @inline function initial_condition_tidal_surge(x, t, equations::ShallowWaterEquations2D)
#     x1, x2 = x
#     b = flat_bathymetry(x1, x2)

#     h = max(equations.threshold_limiter, equations.H0 - b)
#     v1 = 0.0
#     v2 = 0.0
#     return SVector(h, h*v1, h*v2, b)
# end

# @inline function boundary_condition_wave_maker(u_inner, normal_direction::AbstractVector,
#                                                x, t, surface_flux_functions,
#                                                equations::ShallowWaterEquations2D)
#     # Extract the numerical flux functions to compute the conservative and nonconservative
#     # pieces of the approximation
#     surface_flux_function, nonconservative_flux_function = surface_flux_functions

#     # Compute the water height from the wave maker input file data
#     # and then clip to avoid negative water heights and division by zero
#     h_ext = max(equations.threshold_limiter, H_from_wave_maker(t) - u_inner[4])

#     h0 = 25.0
#     v1_ext = 2 * (sqrt(equations.gravity * h_ext) - sqrt(equations.gravity * h0))

#     # Create the external solution state in the conservative variables
#     u_outer = SVector(h_ext, h_ext * v1_ext, zero(eltype(x)), u_inner[4])

#     # Calculate the boundary flux and nonconservative contributions
#     flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

#     noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction,
#                                                  equations)

#     return flux, noncons_flux
# end;






# # export initial_condition_tidal_surge, boundary_condition_tidal_open_left, boundary_condition_tidal_open_right
# # # ---------------------------
# # # Tidal forcing parameters
# # # ---------------------------
# # const T_tide = 12.42 * 3600.0
# # const ω      = 2π / T_tide

# # # Left/right amplitudes and phases
# # const A_left   = 1.5         # [m]
# # const A_right  = 1.2         # [m]  (example: slight damping)
# # const ϕ_left   = 0.0         # [rad]
# # const ϕ_right  = -0.20       # [rad] (example: right lags left)

# # # Optional smooth startup (highly recommended)
# # const T_ramp = 2.0 * 3600.0  # 2 hours
# # @inline ramp(t) = min(1.0, max(0.0, t / T_ramp))

# # @inline tidal_signal(A, ϕ, t) = ramp(t) * A * sin(ω * t + ϕ)

# # @inline eta_left(t, equations::ShallowWaterEquations2D)  = equations.H0 + tidal_signal(A_left,  ϕ_left,  t)
# # @inline eta_right(t, equations::ShallowWaterEquations2D) = equations.H0 + tidal_signal(A_right, ϕ_right, t)


# # # ---------------------------
# # # Initial condition
# # # ---------------------------
# # @inline function initial_condition_tidal_surge(x, t, equations::ShallowWaterEquations2D)
# #     x1, x2 = x
# #     b = flat_bathymetry(x1, x2)

# #     h = max(equations.threshold_limiter, equations.H0 - b)
# #     v1 = 0.0
# #     v2 = 0.0
# #     return SVector(h, h*v1, h*v2, b)
# # end


# # # ---------------------------
# # # Flather ghost-state helper
# # # ---------------------------
# # @inline function flather_outer_state(u_inner, normal_direction, eta_ext, vn_ext,
# #                                      equations::ShallowWaterEquations2D)
# #     h_i, hv1_i, hv2_i, b_i = u_inner
# #     nx, ny = normal_direction

# #     h_floor = max(equations.threshold_limiter, 1e-6)
# #     h_safe  = max(h_i, h_floor)

# #     # interior velocities
# #     v1_i = hv1_i / h_safe
# #     v2_i = hv2_i / h_safe

# #     # interior elevation
# #     eta_i = h_i + b_i

# #     # exterior depth from prescribed elevation
# #     h_o = max(h_floor, eta_ext - b_i)

# #     # wave celerity factor
# #     H_eff = max(0.5 * (h_safe + h_o), h_floor)
# #     cfac  = sqrt(equations.gravity / H_eff)   # sqrt(g/H)

# #     # normal/tangential decomposition
# #     vn_i = v1_i * nx + v2_i * ny
# #     tx, ty = -ny, nx
# #     vt_i = v1_i * tx + v2_i * ty

# #     # Flather relation (ROMS form): vn = vn_ext + sqrt(g/H)*(eta_i - eta_ext)
# #     vn_o = vn_ext + cfac * (eta_i - eta_ext)

# #     # (Optional safety caps while tuning)
# #     vn_o = clamp(vn_o, -10.0, 10.0)
# #     vt_i = clamp(vt_i, -10.0, 10.0)

# #     # reconstruct exterior velocity
# #     v1_o = vn_o * nx + vt_i * tx
# #     v2_o = vn_o * ny + vt_i * ty

# #     return SVector(h_o, h_o * v1_o, h_o * v2_o, b_i)
# # end


# # # ---------------------------
# # # Open boundary conditions
# # # ---------------------------
# # @inline function boundary_condition_tidal_open_left(u_inner, normal_direction,
# #                                                     x, t, surface_flux_functions,
# #                                                     equations::ShallowWaterEquations2D)
# #     sf, ncf = surface_flux_functions

# #     η_ext = eta_left(t, equations)
# #     vn_ext = 0.0  # external barotropic normal velocity (unknown => 0 default)

# #     u_outer = flather_outer_state(u_inner, normal_direction, η_ext, vn_ext, equations)

# #     flux = sf(u_inner, u_outer, normal_direction, equations)
# #     noncons_flux = ncf(u_inner, u_outer, normal_direction, equations)
# #     return flux, noncons_flux
# # end

# # @inline function boundary_condition_tidal_open_right(u_inner, normal_direction,
# #                                                      x, t, surface_flux_functions,
# #                                                      equations::ShallowWaterEquations2D)
# #     sf, ncf = surface_flux_functions

# #     η_ext = eta_right(t, equations)
# #     vn_ext = 0.0

# #     u_outer = flather_outer_state(u_inner, normal_direction, η_ext, vn_ext, equations)

# #     flux = sf(u_inner, u_outer, normal_direction, equations)
# #     noncons_flux = ncf(u_inner, u_outer, normal_direction, equations)
# #     return flux, noncons_flux
# # end

# # export initial_condition_tidal_surge, boundary_condition_tidal_inflow, boundary_condition_tidal_outflow

# # const T_tide = 12.42 * 3600.0
# # const ω      = 2π / T_tide
# # const A      = 1.5  # tidal amplitude [m]

# # @inline tidal_elevation(t) = A * sin(ω * t)

# # @inline function initial_condition_tidal_surge(x, t, equations::ShallowWaterEquations2D)
# #     x1, x2 = x
# #     b = flat_bathymetry(x1, x2)

# #     # still water at eta = H0
# #     h = max(equations.threshold_limiter, equations.H0 - b)
# #     v1 = 0.0
# #     v2 = 0.0
# #     return SVector(h, h * v1, h * v2, b)
# # end

# # # Build ghost/outer state using Flather relation for normal velocity:
# # # vn = vn_ext + sqrt(g/H)*(eta_inner - eta_ext)
# # @inline function flather_outer_state(u_inner, normal_direction, eta_ext, vn_ext,
# #                                      equations::ShallowWaterEquations2D)
# #     h_inner, hv1_inner, hv2_inner, b_inner = u_inner
# #     nx, ny = normal_direction

# #     h_floor = max(equations.threshold_limiter, 1e-6)
# #     h_i = max(h_inner, h_floor)

# #     # interior velocities
# #     v1_i = hv1_inner / h_i
# #     v2_i = hv2_inner / h_i

# #     # interior free-surface elevation
# #     eta_i = h_inner + b_inner

# #     # outer depth from prescribed external elevation
# #     h_o = max(h_floor, eta_ext - b_inner)

# #     # effective depth H in Flather wave speed factor
# #     H_eff = max(0.5 * (h_i + h_o), h_floor)
# #     cfac  = sqrt(equations.gravity / H_eff)  # sqrt(g/H)

# #     # normal/tangential decomposition (outward-normal convention)
# #     vn_i = v1_i * nx + v2_i * ny
# #     tx, ty = -ny, nx
# #     vt_i = v1_i * tx + v2_i * ty

# #     # Flather normal velocity
# #     vn_o = vn_ext + cfac * (eta_i - eta_ext)

# #     # optional safety cap while debugging
# #     vn_o = clamp(vn_o, -10.0, 10.0)
# #     vt_i = clamp(vt_i, -10.0, 10.0)

# #     # reconstruct outer velocity
# #     v1_o = vn_o * nx + vt_i * tx
# #     v2_o = vn_o * ny + vt_i * ty

# #     return SVector(h_o, h_o * v1_o, h_o * v2_o, b_inner)
# # end

# # @inline function boundary_condition_tidal_inflow(u_inner, normal_direction,
# #                                                   x, t, surface_flux_functions,
# #                                                   equations::ShallowWaterEquations2D)
# #     surface_flux_function, nonconservative_flux_function = surface_flux_functions

# #     # External sea level forcing (M2 tide) at LEFT/open boundary
# #     eta_ext = equations.H0 + tidal_elevation(t)
# #     vn_ext  = 0.0  # external barotropic normal velocity (set 0 if unknown)

# #     u_outer = flather_outer_state(u_inner, normal_direction, eta_ext, vn_ext, equations)

# #     flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
# #     noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction, equations)
# #     return flux, noncons_flux
# # end

# # @inline function boundary_condition_tidal_outflow(u_inner, normal_direction,
# #                                                    x, t, surface_flux_functions,
# #                                                    equations::ShallowWaterEquations2D)
# #     surface_flux_function, nonconservative_flux_function = surface_flux_functions

# #     # Radiation/open boundary with reference sea level on RIGHT boundary
# #     eta_ext = equations.H0
# #     vn_ext  = 0.0

# #     u_outer = flather_outer_state(u_inner, normal_direction, eta_ext, vn_ext, equations)

# #     flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)
# #     noncons_flux = nonconservative_flux_function(u_inner, u_outer, normal_direction, equations)
# #     return flux, noncons_flux
# # end

