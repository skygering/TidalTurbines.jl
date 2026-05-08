using HOHQMesh, OrdinaryDiffEqSSPRK, Trixi, TrixiShallowWater, Trixi2Vtk, TidalTurbines
using ForwardDiff

# [... your setup code ...]

# ============================================================================
# STEP 1: Modify build_semi to accept dual numbers in uEltype
# ============================================================================
function build_semi(p_list, base; uEltype=Float64)
    # turbines - this now must handle whatever type p_list contains
    turbines = [Turbine(; x0 = p[1], y0 = p[2], u_rated, u_in, u_out, h_min) for p in p_list]

    # source terms
    sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
    friction_source = source_terms_manning_friction
    turbine_source = make_turbine_source(turbines)
    all_source_terms = combine_source_terms([sponge_source, friction_source, turbine_source])

    # KEY: Pass uEltype to support dual numbers
    semi = SemidiscretizationHyperbolic(
        base.mesh, base.equations, base.initial_condition, base.solver;
        boundary_conditions = base.boundary_condition,
        source_terms = all_source_terms,
        uEltype = uEltype  # ← CRITICAL for AD
    )
    return semi, turbines
end

# ============================================================================
# STEP 2: Helper functions to convert between flat vector and tuple format
# ============================================================================

"""
    flatten_turbine_locations(p_list::Vector{Tuple})

Convert [(x1,y1), (x2,y2), ...] to [x1, y1, x2, y2, ...]
"""
function flatten_turbine_locations(p_list)
    return [p[i] for p in p_list for i in 1:2]
end

"""
    unflatten_turbine_locations(p_flat::Vector, n_turbines::Int)

Convert [x1, y1, x2, y2, ...] back to [(x1,y1), (x2,y2), ...]
"""
function unflatten_turbine_locations(p_flat, n_turbines)
    return [Tuple(p_flat[2i-1:2i]) for i in 1:n_turbines]
end



# ============================================================================
# STEP 3: Create an objective function that works with flat vectors and AD
# ============================================================================

function objective_ad(p_flat::AbstractVector, base, n_turbines; 
                      tspan, saveat=0.05, rho=1000.0)
    
    # Unflatten parameters back to list of locations
    p_list = unflatten_turbine_locations(p_flat, n_turbines)
    
    # Get element type (handles both Float64 and ForwardDiff.Dual)
    element_type = typeof(first(first(p_list)))
    
    # Build semidiscretization with appropriate uEltype
    semi, turbines = build_semi(p_list, base; uEltype=element_type)
    ode = semidiscretize(semi, tspan)

    callbacks = CallbackSet()
    stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

    sol = solve(ode, SSPRK43(stage_limiter!);
                dt = Δt,
                ode_default_options()...,
                callback = callbacks,
                adaptive = false,
                saveat)

    E, _ = compute_total_turbine_energy(sol, semi, turbines; rho)
    return -E
end

# ============================================================================
# STEP 4: Compute the gradient
# ============================================================================

base = build_base_simulation()

# Define initial turbine locations
p_initial = [(Lx / 2, 0.25 * Ly)]  # List of tuples
n_turbines = length(p_initial)

# Flatten for ForwardDiff
p_flat = flatten_turbine_locations(p_initial)

## Create a closure with fixed arguments
objective_closure(p) = objective_ad(p, base, n_turbines; tspan=(0.0, 2.0))

# Compute objective value at initial point
J0 = objective_closure(p_flat)
println("Objective value: $J0")

# Compute gradient
∇J = ForwardDiff.gradient(objective_closure, p_flat)
println("Gradient w.r.t. turbine locations: $∇J")

# Unflatten gradient for interpretation
∇J_list = unflatten_turbine_locations(∇J, n_turbines)
println("Gradient per turbine (∂J/∂x, ∂J/∂y):")
for (i, grad) in enumerate(∇J_list)
    println("  Turbine $i: ∇J = $grad")
end

## ============================================================================
# STEP 5: (Optional) Use for simple optimization
# ============================================================================

using Optim

result = optimize(objective_closure, p_flat, BFGS())

p_optimal = Optim.minimizer(result)
p_optimal_list = unflatten_turbine_locations(p_optimal, n_turbines)

println("Optimal turbine locations: $p_optimal_list")
println("Optimal objective value: $(Optim.minimum(result))")

## Testing ground
using ForwardDiff

# Simple dual number test
x_dual = ForwardDiff.Dual(2.5, 1.0)
y_dual = ForwardDiff.Dual(3.0, 0.0)

# Try step by step
println("x_dual: $x_dual, type: $(typeof(x_dual))")
println("y_dual: $y_dual, type: $(typeof(y_dual))")

# Test promote_type
T = promote_type(typeof(x_dual), typeof(y_dual))
println("Promoted type: $T")

# Try conversion
try
    x_converted = T(x_dual)
    println("✓ x_converted: $x_converted")
catch e
    println("✗ Conversion failed: $e")
end

## Debug details to find nan
using ForwardDiff

function objective_ad_debug_detailed(p_flat::AbstractVector, base, n_turbines; 
                      tspan, saveat=0.05, rho=1000.0)
    
    println("=== objective_ad_debug_detailed ===")
    println("Input p_flat: $p_flat")
    
    p_list = unflatten_turbine_locations(p_flat, n_turbines)
    println("✓ Unflattened p_list: $p_list")
    
    element_type = typeof(first(first(p_list)))
    println("Element type: $element_type")
    
    semi, turbines = build_semi(p_list, base; uEltype=element_type)
    println("✓ Built semi")
    
    ode = semidiscretize(semi, tspan)
    println("✓ Semidiscretized")

    callbacks = CallbackSet()
    stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

    println("About to solve ODE...")
    sol = solve(ode, SSPRK43(stage_limiter!);
                dt = Δt,
                ode_default_options()...,
                callback = callbacks,
                adaptive = false,
                saveat)


    # After the ODE solve in the debug code:
    println("\n=== Diagnosing compute_total_turbine_energy ===")

    # Try calling it step by step
    sol_end = sol.u[end]
    println("sol_end type: $(typeof(sol_end))")
    println("First element: $(sol_end[1])")

    # Try computing energy with Float64 first
    E_float, _ = compute_total_turbine_energy(sol, semi, turbines; rho)
    println("Energy (Float64): $E_float")

    # Now the issue - try getting individual turbine contributions
    for (i, turb) in enumerate(turbines)
        println("\nTurbine $i at ($(turb.x0), $(turb.y0))")
        # Try to manually compute power at a point
        try
            x_test = SVector(turb.x0, turb.y0)
            power_dens = turbine_power_density_single(sol_end, x_test, semi.equations, turb; rho)
            println("  Power density: $power_dens, is NaN in partials? $(isnan(power_dens.partials[1]))")
        catch e
            println("  Error: $e")
        end
    end            
    
    println("✓ ODE solved")
    println("sol.u[end] type: $(typeof(sol.u[end]))")
    println("sol.u[end] contains NaNs? $(any(isnan.(sol.u[end])))")

    E, _ = compute_total_turbine_energy(sol, semi, turbines; rho)
    println("Energy computed: $E")
    println("Energy is NaN? $(isnan(E))")
    
    result = -E
    println("Final result: $result, Is NaN? $(isnan(result))")
    
    return result
end

# Test 1: Works with Float64?
println("=== TEST 1: Regular Float64 ===")
p_flat_float = flatten_turbine_locations([(Lx / 2, 0.25 * Ly)])
objective_closure_debug(p) = objective_ad_debug_detailed(p, base, n_turbines; tspan=(0.0, 2.0))
J0 = objective_closure_debug(p_flat_float)

# Test 2: Try with single dual
println("\n=== TEST 2: Single Dual Number ===")
p_flat_dual = [ForwardDiff.Dual(Lx/2, 1.0), ForwardDiff.Dual(0.25*Ly, 0.0)]
try
    J_dual = objective_closure_debug(p_flat_dual)
    println("Result with dual: $J_dual")
catch e
    println("ERROR: $e")
    import Traceback
    showerror(stderr, e)
end

##
println("\n=== TEST 4: Minimal ODE solve with duals ===")

p_flat_dual = [ForwardDiff.Dual(Lx/2, 1.0), ForwardDiff.Dual(0.25*Ly, 0.0)]
p_list = unflatten_turbine_locations(p_flat_dual, 1)
element_type = typeof(first(first(p_list)))

# Build semi WITHOUT any source terms
mesh = UnstructuredMesh2D(joinpath("examples/tidal_headland_objective", "tidal_headland.mesh"))
basis = LobattoLegendreBasis(3)
surface_flux = (
    FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
    flux_nonconservative_chen_noelle,
)
volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)
equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)

bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)
initial_condition = make_initial_condition_tidal_surge(bathy)

wm_left  = WaveMaker(Aₜ, ω)
wm_right = WaveMaker(-Aₜ, ω)
bc_left  = BoundaryConditionDirichlet(make_dirichlet_state(wm_left, bathy))
bc_right = BoundaryConditionDirichlet(make_dirichlet_state(wm_right, bathy))
boundary_condition = (;
    Bottom = boundary_condition_slip_wall,
    Top    = boundary_condition_slip_wall,
    Left   = bc_left,
    Right  = bc_right,
)

# NO SOURCE TERMS
semi_minimal = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    uEltype = element_type
)

ode = semidiscretize(semi_minimal, (0.0, 0.5))

println("Solving minimal ODE without source terms...")
try
    sol = solve(ode, SSPRK43();
                dt = Δt,
                ode_default_options()...,
                adaptive = false,
                saveat = 0.05)
    
    println("✓ Solve succeeded")
    println("sol.u[end][1]: $(sol.u[end][1])")
    if typeof(sol.u[end][1]) <: ForwardDiff.Dual
        println("Has NaN? $(isnan(sol.u[end][1].partials[1]))")
    end
catch e
    println("✗ Error: $e")
    showerror(stderr, e)
end

##
println("\n=== TEST 5a: Add sponge source only ===")

sponge_source = make_sponge_source(; Lx, σ_max=σₘ)

semi_sponge = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = sponge_source,
    uEltype = element_type
)

ode = semidiscretize(semi_sponge, (0.0, 0.5))

println("Solving with sponge source only...")
try
    sol = solve(ode, SSPRK43();
                dt = Δt,
                ode_default_options()...,
                adaptive = false,
                saveat = 0.05)
    
    println("✓ Sponge source OK")
    if typeof(sol.u[end][1]) <: ForwardDiff.Dual
        println("Has NaN? $(isnan(sol.u[end][1].partials[1]))")
    end
catch e
    println("✗ Sponge source ERROR: $e")
end

println("\n=== TEST 5b: Add friction source only ===")

friction_source = source_terms_manning_friction

semi_friction = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = friction_source,
    uEltype = element_type
)

ode = semidiscretize(semi_friction, (0.0, 0.5))

println("Solving with friction source only...")
try
    sol = solve(ode, SSPRK43();
                dt = Δt,
                ode_default_options()...,
                adaptive = false,
                saveat = 0.05)
    
    println("✓ Friction source OK")
    if typeof(sol.u[end][1]) <: ForwardDiff.Dual
        println("Has NaN? $(isnan(sol.u[end][1].partials[1]))")
    end
catch e
    println("✗ Friction source ERROR: $e")
end

println("\n=== TEST 5c: Add turbine source only ===")

turbines = [Turbine(; x0 = p[1], y0 = p[2], u_rated, u_in, u_out, h_min) for p in p_list]
turbine_source = make_turbine_source(turbines)

semi_turbine = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = turbine_source,
    uEltype = element_type
)

ode = semidiscretize(semi_turbine, (0.0, 0.5))

println("Solving with turbine source only...")
try
    sol = solve(ode, SSPRK43();
                dt = Δt,
                ode_default_options()...,
                adaptive = false,
                saveat = 0.05)
    
    println("✓ Turbine source OK")
    if typeof(sol.u[end][1]) <: ForwardDiff.Dual
        println("Has NaN? $(isnan(sol.u[end][1].partials[1]))")
    end
catch e
    println("✗ Turbine source ERROR: $e")
end

##
using ForwardDiff

# Test the friction function directly with a dual
h_dual = ForwardDiff.Dual(0.1, 0.01)  # Small h value
hv_1_dual = ForwardDiff.Dual(0.5, 0.0)
hv_2_dual = ForwardDiff.Dual(0.3, 0.0)

u_test = SVector(h_dual, hv_1_dual, hv_2_dual, ForwardDiff.Dual(0.0, 0.0))

equations_test = ShallowWaterEquations2D(gravity=1.0, H0=0.0)

result = source_terms_manning_friction(u_test, SVector(0.0, 0.0), 0.0, equations_test)

println("Result: $result")
println("Result[2]: $(result[2])")
if typeof(result[2]) <: ForwardDiff.Dual
    println("Partial: $(result[2].partials)")
    println("Has NaN? $(isnan(result[2].partials[1]))")
end

##
using Statistics
println("\n=== TEST 5b: Add friction source only - CHECK MIN H ===")

p_flat_dual = [ForwardDiff.Dual(Lx/2, 1.0), ForwardDiff.Dual(0.25*Ly, 0.0)]
p_list = unflatten_turbine_locations(p_flat_dual, 1)
element_type = typeof(first(first(p_list)))

friction_source = source_terms_manning_friction

semi_friction = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = friction_source,
    uEltype = element_type
)

ode = semidiscretize(semi_friction, (0.0, 0.5))

println("Solving with friction source...")
sol = solve(ode, SSPRK43();
            dt = Δt,
            ode_default_options()...,
            adaptive = false,
            saveat = 0.05)

println("\n=== Analyzing h values throughout solution ===")
for (k, u) in enumerate(sol.u)
    # Extract h values (every 4th element starting at 1)
    h_values = u[1:4:end]
    
    # Get real parts (strip dual numbers)
    h_real = real.(h_values)
    
    min_h = minimum(h_real)
    max_h = maximum(h_real)
    mean_h = mean(h_real)
    
    println("Step $k (t=$(sol.t[k])): h ∈ [$min_h, $max_h], mean=$mean_h")
    
    # Check for NaN - simpler approach
    has_nan = any(isnan.(h_values))
    if has_nan
        println("  ⚠️ NaN in h values detected!")
        
        # Find which ones are NaN
        for (i, h_val) in enumerate(h_values)
            if isnan(h_val)
                println("    Node $i: $h_val")
            end
        end
    end
end

# Calculate global minimum h
all_h_real = reduce(vcat, [real.(u[1:4:end]) for u in sol.u])
println("\nGlobal minimum h: $(minimum(all_h_real))")
println("Global maximum h: $(maximum(all_h_real))")

##
println("\n=== TEST 5d: Friction + Sponge combined ===")

p_flat_dual = [ForwardDiff.Dual(Lx/2, 1.0), ForwardDiff.Dual(0.25*Ly, 0.0)]
p_list = unflatten_turbine_locations(p_flat_dual, 1)
element_type = typeof(first(first(p_list)))

sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
friction_source = source_terms_manning_friction
all_source_terms = combine_source_terms([sponge_source, friction_source])

semi_combined = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = all_source_terms,
    uEltype = element_type
)

ode = semidiscretize(semi_combined, (0.0, 0.5))

println("Solving with sponge + friction...")
sol = solve(ode, SSPRK43();
            dt = Δt,
            ode_default_options()...,
            adaptive = false,
            saveat = 0.05)

println("\nStep 1 (t=0.0): h[1] = $(sol.u[1][1])")
println("Step 2 (t=0.05): h[1] = $(sol.u[2][1])")
println("Has NaN in step 2? $(isnan(sol.u[2][1]))")

##
println("\n=== TEST 5e: Sponge ONLY ===")

friction_source = source_terms_manning_friction


semi_sponge_only = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver;
    boundary_conditions = boundary_condition,
    source_terms = friction_source,
    uEltype = element_type
)

ode = semidiscretize(semi_sponge_only, (0.0, 0.5))

println("Solving with sponge only...")
sol = solve(ode, RK4();
            dt = Δt,
            ode_default_options()...,
            adaptive = false,
            saveat = 0.05)
println("Semi cache type: $(typeof(semi_sponge_only.cache))")
println("Semi has limiter? $(hasproperty(semi_sponge_only, :limiter))")
println("Step 2: h[1] = $(sol.u[2][1])")
println("Has NaN? $(isnan(sol.u[2][1]))")