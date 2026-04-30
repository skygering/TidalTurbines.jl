# Optimization framework

# Load packages
using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk
using ForwardDiff

## -----------------------------------------
# Use shallow water equations
#  -----------------------------------------
equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)


## -----------------------------------------
# Initial turbine position and radius
#  -----------------------------------------
locs = [11.0, 4.0]#, 3.0, 2.0, 2.0, 3.0];  # (x1, y1, x2, y2)
Diameter = 1.0;
R = 0.5;

## ----------------------------------------- 
# Bathymetry & Mesh
#  -----------------------------------------
#flat_bathymetry(x, y) = -2.0
include("bathymetry.jl")
include("boundary_conditions.jl")
bathy = FlatBathymetry(-2.0)


## Create an unstructured mesh project
flat_bottom = newProject("flat_bottom", "out");

setPolynomialOrder!(flat_bottom, 1)
setMeshFileFormat!(flat_bottom, "ISM-V2");
HOHQMesh.getModelDict(flat_bottom);

# Set domain size (edges) we can set a background Cartesian box mesh required to define order `[top, left, bottom, right]`.
bounds = [8.0, 0.0, 0.0, 22.0]
N = [22, 8, 0] # Coarse background grid: 8 x 4 in (x, y)
addBackgroundGrid!(flat_bottom, bounds, N)

# Add refinement regions around turbines
for i in 1:2:length(locs)
    name = "turbine_$(i÷2 + 1)"
    wake_refinement = newRefinementLine(name, "smooth", [locs[i]-2*R, locs[i+1], 0.0], 
                                        [locs[i]+5*R, locs[i+1], 0.0], 0.5*R, 3*R);
    add!(flat_bottom, wake_refinement)
end

# Plot and save mesh
plotProject!(flat_bottom, GRID + REFINEMENTS);

generate_mesh(flat_bottom);

mesh_file = joinpath(@__DIR__, "..", "out", "flat_bottom.mesh")
mesh = UnstructuredMesh2D(mesh_file)


## ----------------------------------------- 
# Initial Conditions
#  -----------------------------------------
initial_condition = make_initial_condition_tidal_surge(bathy)
@inline function initial_condition_tsunami(x, t, equations::ShallowWaterEquations2D)
    ## Initially water is at rest
    v1 = 0.0
    v2 = 0.0

    ## Bottom topography values are computed from the bicubic spline created above
    x1, x2 = x
    b = flat_bathymetry(x1, x2)

    h = max(equations.threshold_limiter, equations.H0 - b) # ensure positive wave height

    ## Return the conservative variables
    return SVector(h, h * v1, h * v2, b)
end

# initial_condition = initial_condition_tsunami;

## ----------------------------------------- 
# Boundary Conditions
#  -----------------------------------------
# RightBC = BoundaryConditionWaterHeight(t -> h_boundary(t), equations)
include("boundary_conditions.jl")

# boundary_condition = (; Bottom = boundary_condition_slip_wall,
#                       Top = boundary_condition_slip_wall,
#                       Right = boundary_condition_slip_wall,
#                       Left = boundary_condition_tidal_surge);
wm = WaveMaker(0.5, 0.25) # amplitude, frequency 
equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)
initial_condition = make_initial_condition_tidal_surge(bathy)
boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = make_boundary_condition_flather_left(wm, bathy),
    Right  = make_boundary_condition_flather_right(bathy),
)

## ----------------------------------------- 
# Turbine forcing & total drag
#  -----------------------------------------

# Define turbine bump function for Ct
@inline function bump_1d(ξ, r)
    """
    1D smooth bump function, sharp and localized.
    Centered at ξ=0 with characteristic radius r.
    """
    s = abs(ξ) / r
    return exp(-s^2 / 0.5)  # tighter Gaussian
end

@inline function bump_2d(x, y, x_turb, y_turb, r)
    """2D separable bump function: ψ(x) * ψ(y)"""
    return bump_1d(x - x_turb, r) * bump_1d(y - y_turb, r)
end

# Define bottom friction with Manning friction + turbine bump(s): Function of turbine locations
@inline function source_terms_all(u, x, t, equations::ShallowWaterEquations2D, locs)

    h, hv_1, hv_2, _ = u
    x_pos, y_pos = x
    
    K_strength = 0.01;        # [K₁, K₂, ...]
    r_turbine = 0.3;       # bump radius

    T = eltype(x)
    
    # ===== Manning Friction =====
    n = 0.001  # friction coefficient
    h_eff = (h^2 + max(h^2, 1e-8)) / (2 * h)  # desingularization
    Sf = -equations.gravity * n^2 * h_eff^(-7/3) * sqrt(hv_1^2 + hv_2^2)
    
    S_hv1 = Sf * hv_1
    S_hv2 = Sf * hv_2
    
    # ===== Turbine Drag =====
    U = sqrt(hv_1^2 + hv_2^2) / max(h_eff, 1e-8)  # flow speed
    n_turbines = div(length(locs), 2)
    
    for i in 1:n_turbines
        x_i = locs[2*i - 1]
        y_i = locs[2*i]
        K_i = K_strength
        
        # Spatial drag field: K_i(x,y) = K_i * ψ_{x_i,r}(x) * ψ_{y_i,r}(y)
        K_field = K_i * bump_2d(x_pos, y_pos, x_i, y_i, r_turbine)
        
        # Drag force: -K_field * U * v
        drag_coeff = K_field * U
        S_hv1 -= drag_coeff * hv_1
        S_hv2 -= drag_coeff * hv_2
    end
    
    return SVector(zero(T), S_hv1, S_hv2, zero(T))
end

# Alternatively, use just bottom friction
@inline function source_terms_manning_friction(u, x, t,
                                               equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001 # friction coefficient
    h = (h^2 + max(h^2, 1e-8)) / (2 * h) # desingularization procedure

    ## Compute the common friction term
    Sf = -equations.gravity * n^2 * h^(-7 / 3) * sqrt(hv_1^2 + hv_2^2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end;
## ----------------------------------------- 
# Weixuan's source terms
#  -----------------------------------------
@inline function Cp_from_Ct(Ct)
    Ct_clamped = clamp(Ct, 0.0, 1.0)
    return 0.5 * (1.0 + sqrt(1.0 - Ct_clamped)) * Ct_clamped
end


const XI_BUMP = 1.45661
const Ct_rated_val = 0.516
const Cp_rated_val = Cp_from_Ct(Ct_rated_val)

const turbines = [
    (
        x0 = 11.0,
        y0 = 4.0,
        D = 1.0,
        r_support = 0.5,
        Ct_rated = Ct_rated_val,
        Cp_rated = Cp_rated_val,
        u_in = 1/14,
        u_rated = 3/14,
        u_out = 6/14,
        h_min = 0.00
    )
]
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

@inline function Ct_turbine(U, turb)
    u_in      = turb.u_in       # 1.0/14
    u_rated   = turb.u_rated    # 3.05/14
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


## ----------------------------------------- 
# Finite volume solver
#  -----------------------------------------
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(4) # polynomial approximation space with degree 7

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight)
# volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
#                                                  volume_flux_dg = volume_flux,
#                                                  volume_flux_fv = surface_flux)

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

## ----------------------------------------- 
# Semidiscretization object
#  -----------------------------------------
#source_terms_inner = (u, x, t, equations) -> source_terms_all(u, x, t, equations, locs);
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_terms_manning_plus_turbines);

tspan = (0.0, 30.0)
ode = semidiscretize(semi, tspan);

## ----------------------------------------- 
# Callbacks
#  -----------------------------------------
analysis_interval = 1000 
analysis_callback = AnalysisCallback(semi, interval = analysis_interval);
save_solution = SaveSolutionCallback(dt = 0.1,
                                     save_initial_solution = true,
                                     save_final_solution = true);
stepsize_callback = StepsizeCallback(cfl = 0.6);
callbacks = CallbackSet(analysis_callback,
                        stepsize_callback,
                        save_solution);
# callbacks = CallbackSet();

## ----------------------------------------- 
# Try running the problem
#  -----------------------------------------
stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))
sol = solve(ode, SSPRK43(stage_limiter!); dt = 0.0001,
            ode_default_options()..., callback = callbacks, adaptive = false);
# ----------------------------------------- 
# Save the solution as a vtk
#  -----------------------------------------
trixi2vtk("out/solution_*.h5", output_directory = "out")


## ----------------------------------------- 
# AD Attempt
#  -----------------------------------------
# doesn't allow callbacks right now...

# jacobian_ad_forward(semi)
u0_ode = Trixi.compute_coefficients(0.0, semi)


J = ForwardDiff.jacobian((du_ode, locs) -> begin
                            equations_inner = equations
                            source_terms_inner = source_terms_manning_plus_turbines

                            semi_inner = Trixi.remake(semi, 
                                source_terms = source_terms_inner,
                                equations = equations_inner, 
                                uEltype = eltype(locs)
                                )
                             Trixi.rhs!(du_ode, u0_ode, semi_inner, 0.0)
                         end, similar(u0_ode), [11.0, 4.0]); 

