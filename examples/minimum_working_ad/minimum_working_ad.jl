using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using ForwardDiff

# ============================================================================
# MINIMAL SHALLOW WATER AD TEST
# ============================================================================

output_directory = "examples/minimum_working_ad/out"
rm(output_directory, recursive=true, force=true)
mkpath(output_directory)

# Simple domain setup
# Lx = 10.0
# Ly = 5.0
# g = 9.81

# Create a simple mesh
# ad_test = newProject("min_ad_test", "examples/minimum_working_ad")
# setPolynomialOrder!(ad_test, 1)
# setMeshFileFormat!(ad_test, "ISM-V2")

# bounds = [Ly, 0.0, 0.0, Lx]
# N = [20, 10, 0]
# addBackgroundGrid!(ad_test, bounds, N)
# generate_mesh(ad_test)

# Plot and save mesh
# plotProject!(ad_test, GRID + REFINEMENTS);
# generate_mesh(ad_test);

# mesh_file = joinpath(@__DIR__, "min_ad_test.mesh")
# mesh = UnstructuredMesh2D(mesh_file)

# bathy = FlatBathymetry(-2.0)

# Create an unstructured mesh project
refined_mesh = newProject("refined_mesh", "examples/minimum_working_ad")
setPolynomialOrder!(refined_mesh, 1)
setMeshFileFormat!(refined_mesh, "ISM-V2");
HOHQMesh.getModelDict(refined_mesh);

# Set domain size (edges) we can set a background Cartesian box mesh required to define order `[top, left, bottom, right]`.
bounds = [4.0, 0.0, 0.0, 10.0]
N = [60, 20, 0] # Coarse background grid: 8 x 4 in (x, y)
addBackgroundGrid!(refined_mesh, bounds, N)

p0 = [3.0, 2.0]
R = 0.25

# Add refinement regions around turbines
for i in 1:2:length(p0)
    name = "turbine_$(i÷2 + 1)"
    wake_refinement = newRefinementLine(name, "smooth", [p0[i]-2*R, p0[i+1], 0.0], 
                                        [p0[i]+5*R, p0[i+1], 0.0], 0.5*R, 3*R);
    add!(refined_mesh, wake_refinement)
end

# Plot and save mesh
plotProject!(refined_mesh, GRID + REFINEMENTS);
generate_mesh(refined_mesh);

mesh_file = joinpath(@__DIR__, "refined_mesh.mesh")
mesh = UnstructuredMesh2D(mesh_file)

## ============================================================================
# Manning friction source term
# ============================================================================

@inline function source_terms_manning_friction(u, x, t, equations::ShallowWaterEquations2D)
    h, hv_1, hv_2, _ = u

    n = 0.001
    h = (h^2 + max(h^2, 1e-8)) / (2 * h)

    Sf = -equations.gravity * n^2 * h^(-7 / 3) * (hv_1^2 + hv_2^2)^(1/2)

    return SVector(zero(eltype(x)), Sf * hv_1, Sf * hv_2, zero(eltype(x)))
end


## ====This worked (didn't yield NaNs in gradient)============================
# Build base simulation
# ============================================================================

function build_shallow_water_simulation()
    mesh_file = joinpath("examples/minimum_working_ad", "refined_mesh.mesh")
    mesh = UnstructuredMesh2D(mesh_file)
    
    basis = LobattoLegendreBasis(4)
    surface_flux = (
        FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
        flux_nonconservative_chen_noelle,
    )
    volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)
    
    equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)
    
    # Simple initial condition: Gaussian bump
    function initial_condition_bump(x, t, equations)
        x_c, y_c = Lx/2, Ly/2
        h = 1.0 + 0.2 * exp(-((x[1] - x_c)^2 + (x[2] - y_c)^2) / 0.5)
        return SVector(h, 0.0, 0.0, 0.0)
    end
    
    # Use slip wall for unstructured mesh
    bc_left  = BoundaryConditionDirichlet(make_dirichlet_state(wm_left, bathy))
    bc_right = BoundaryConditionDirichlet(make_dirichlet_state(wm_right, bathy))
    boundary_condition = (;
        Bottom = boundary_condition_slip_wall,
        Top    = boundary_condition_slip_wall,
        Left   = bc_left,
        Right  = bc_right,
    )
    
    return (;
        mesh,
        solver,
        equations,
        initial_condition = initial_condition_bump,
        boundary_condition,
    )
end

# ============================================================================
# Simple objective: total water volume at final time
# ============================================================================

function water_volume(u_ode, semi)
    """Integrate total water height (volume)"""
    u = Trixi.wrap_array(u_ode, semi)
    solver = semi.solver
    cache = semi.cache
    
    volume = zero(eltype(u))
    
    for element in eachelement(solver, cache)
        for j in eachnode(solver), i in eachnode(solver)
            h = u[1, i, j, element]
            w_i = solver.basis.weights[i]
            w_j = solver.basis.weights[j]
            J = cache.elements.inverse_jacobian[i, j, element]
            dV = w_i * w_j / J
            
            volume += h * dV
        end
    end
    
    return volume
end

# ============================================================================
# Objective function parameterized by peak height
# ============================================================================

function objective_height(peak_height; tspan=(0.0, 1.0))
    """
    Objective: integrate final water volume with respect to initial peak height
    """
    base = build_shallow_water_simulation()
    
    # Modify initial condition to use peak_height parameter
    function initial_condition_param(x, t, equations)
        x_c, y_c = Lx/2, Ly/2
        h = 1.0 + peak_height * exp(-((x[1] - x_c)^2 + (x[2] - y_c)^2) / 0.5)
        return SVector(h, 0.0, 0.0, 0.0)
    end
    
    semi = SemidiscretizationHyperbolic(
        base.mesh, base.equations, initial_condition_param, base.solver;
        boundary_conditions = base.boundary_condition,
        source_terms = source_terms_manning_friction,
        uEltype = typeof(peak_height)  # CRITICAL: Enable duals
    )
    
    ode = semidiscretize(semi, tspan)
    
    # Simple RK4 solver (no stage limiter, more AD-friendly)
    sol = solve(ode, RK4();
                dt = 0.001,
                adaptive = false,
                save_everystep = false)
    
    # Return final water volume
    return water_volume(sol.u[end], semi)
end

# ============================================================================
# Test 1: Objective works with Float64
# ============================================================================
tspan = (0,0.5)
println("=== TEST 1: Objective with Float64 ===")
peak_height_0 = 0.2
V_0 = objective_height(peak_height_0; tspan=tspan)
println("Water volume at peak_height=$peak_height_0: $V_0")

# ============================================================================
# Test 2: Try with ForwardDiff.Dual
# ============================================================================

println("\n=== TEST 2: Objective with Dual number ===")
peak_height_dual = ForwardDiff.Dual(0.2, 1.0)
try
    V_dual = objective_height(peak_height_dual; tspan=tspan)
    println("Result: $V_dual")
    println("Partial derivative: $(V_dual.partials[1])")
    
    if isnan(V_dual.partials[1])
        println("⚠️ NaN in partial!")
    else
        println("✓ Valid gradient!")
    end
catch e
    println("✗ Error: $e")
    showerror(stderr, e)
end

# ============================================================================
# Test 3: Try ForwardDiff.derivative
# ============================================================================

println("\n=== TEST 3: ForwardDiff.derivative ===")
try
    dV_dp = ForwardDiff.derivative(objective_height, 0.5)
    println("Derivative: $dV_dp")
    println("✓ Success!")
catch e
    println("✗ Error: $e")
end

## ============================================================================
# Test 2b: SIMPLIFIED - Just check if duals stay healthy
# ============================================================================

println("\n=== TEST 2b: Check dual health at select times ===")

function objective_height_check(peak_height; tspan=(0.0,1.0))
    base = build_shallow_water_simulation()
    
    function initial_condition_param(x, t, equations)
        x_c, y_c = Lx/2, Ly/2
        h = 1.0 + peak_height * exp(-((x[1] - x_c)^2 + (x[2] - y_c)^2) / 0.5)
        return SVector(h, 0.0, 0.0, 0.0)
    end
    
    semi = SemidiscretizationHyperbolic(
        base.mesh, base.equations, initial_condition_param, base.solver;
        boundary_conditions = base.boundary_condition,
        uEltype = typeof(peak_height)
    )
    
    ode = semidiscretize(semi, tspan)
    
    # Save at specific times for inspection
    sol = solve(ode, RK4();
                dt = 0.01,
                adaptive = false,
                saveat = 0.1)  # Save every 0.1 time units
    
    # Check all saved time steps
    for (step, (t, u_ode)) in enumerate(zip(sol.t, sol.u))
        u = Trixi.wrap_array(u_ode, semi)
        
        # Get statistics
        h_vals = u[1, :, :, :]
        hv_1_vals = u[2, :, :, :]
        hv_2_vals = u[3, :, :, :]
        
        h_min = minimum(real.(h_vals))
        h_max = maximum(real.(h_vals))
        hv1_min = minimum(real.(hv_1_vals))
        hv1_max = maximum(real.(hv_1_vals))
        hv2_min = minimum(real.(hv_2_vals))
        hv2_max = maximum(real.(hv_2_vals))
        
        println("Step $step (t=$t):")
        println("  h:   [$h_min, $h_max]")
        println("  hv_1: [$hv1_min, $hv1_max]")
        println("  hv_2: [$hv2_min, $hv2_max]")
        
        # Check for NaN in partials
        h_has_nan = any(isnan.(h_vals))
        hv1_has_nan = any(isnan.(hv_1_vals))
        hv2_has_nan = any(isnan.(hv_2_vals))
        
        if h_has_nan || hv1_has_nan || hv2_has_nan
            println("  ⚠️ NaN DETECTED in partials!")
        else
            println("  ✓ All duals healthy")
        end
    end
    
    return water_volume(sol.u[end], semi)
end

peak_height_dual = ForwardDiff.Dual(0.2, 1.0)
V_dual = objective_height_check(peak_height_dual; tspan=tspan)
println("\nFinal Result: $V_dual")