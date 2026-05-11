using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using Trixi2Vtk
using TidalTurbines
using ForwardDiff

#= 
This setup is based on the paper:

Jordan et al., Combining shallow-water and analytical wake models for tidal array micro-siting, 
Journal of Ocean Engineering and Marine Energy (2022)

https://link.springer.com/article/10.1007/s40722-022-00225-2
=#

output_directory = "examples/tidal_headland_objective/out"
rm(output_directory, recursive=true, force=true)
mkpath(output_directory)

# problem setup + non_dimensionalization
g = 9.81 # m/s^2
D = 20   # m

# domain values
Lx = 1280 / D
Ly = 480 / D
Lₛ = 0.1 * Lx # sponge layer length
H₀ = -50 / D # bathymetry depth
Hₛ = -5 / D # headland/shelf depth
Rₛ = 160 / D # headland/shelf radius

# tidal values
T = 1 * 60 * 60 * sqrt(g / D)
ω = 2π / T
Aₜ = 4 * 0.275 / D

# turbine values
u_in = 1 / sqrt(D * g) # cut-in speed
u_rated = 3 / sqrt(D * g) # approx rated velocity
u_out = 5 / sqrt(D * g) # cut-in speed
h_min = 0.02 / D

# sponge layer strength σ [1/s]
# time scale: t ~ L_sponge / u_rated [s]
# non-dimensional control paaram: tσ ~5 is good absorbing layer
σₘ = 5 * u_rated / Lₛ

# solver time values
t_out = round(Int, T / 10)
Δt = 0.5 * sqrt(g / D) / 4
tspan = (0.0, 0.25 * T)

# (1) create mesh
tidal = newProject("tidal_headland", "examples/tidal_headland_objective")
setPolynomialOrder!(tidal, 1)
setMeshFileFormat!(tidal, "ISM-V2")
HOHQMesh.getModelDict(tidal)

bounds = [Ly, 0.0, 0.0, Lx] # [top, left, bottom, right]
N = [32, 16, 0]
addBackgroundGrid!(tidal, bounds, N)
headland_core = newRefinementCenter(
    "headland_core", "smooth",
    [Lx/2, Ly, 0.0],
    Lx/(32*2),
    1.5Rₛ
)
add!(tidal, headland_core)
channel = newRefinementLine("channel", "smooth", [15.0, 7.5, 0.0],
                         [45.0, 7.5, 0.0], Lx/(32*2), 10.0)
add!(tidal, channel)
generate_mesh(tidal)

# function that returns solver setup
function build_base_simulation()
    # mesh
    mesh_file = joinpath("examples/tidal_headland_objective", "tidal_headland.mesh")
    mesh = UnstructuredMesh2D(mesh_file)
    
    # solver
    basis = LobattoLegendreBasis(3)
    surface_flux = (
        FluxHydrostaticReconstruction(flux_hll_chen_noelle, hydrostatic_reconstruction_chen_noelle),
        flux_nonconservative_chen_noelle,
    )
    volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)
    
    # equations
    equations = ShallowWaterEquations2D(gravity = 1.0, H0 = 0.0)
    
    # intitial conditions
    bathy = HeadlandBathymetry(Lx/2, Ly, Rₛ, H₀, Hₛ)
    initial_condition = make_initial_condition_tidal_surge(bathy)
    
    # boundary conditions
    wm_left  = WaveMaker(Aₜ, ω)
    wm_right = WaveMaker(-Aₜ, ω)
    bc_left  = DirichletWaveBC(wm_left, bathy)
    bc_right = DirichletWaveBC(wm_right, bathy)
    boundary_condition = (;
        Bottom = boundary_condition_slip_wall,
        Top    = boundary_condition_slip_wall,
        Left   = bc_left,
        Right  = bc_right,
    )
    return (; mesh, solver, equations, initial_condition, boundary_condition)
end

function build_semi(p, base)
    # p is a flat vector: [x1, y1, x2, y2, ...] for n turbines
    n = length(p) ÷ 2
    turbines = [Turbine(; x0 = p[2i-1], y0 = p[2i], u_rated, u_in, u_out, h_min) for i in 1:n]

    # turbines
    # turbines = [Turbine(; x0 = p[1], y0 = p[2], u_rated, u_in, u_out, h_min) for p in p_list]

    # source terms
    friction_source = FrictionSource()
    sponge_source = SpongeSource(Lx, σₘ)
    turbine_source = TurbineFriction(turbines)
    source_terms = CombinedSource((sponge_source, friction_source, turbine_source))

    semi = SemidiscretizationHyperbolic(
        base.mesh, base.equations, base.initial_condition, base.solver;
        boundary_conditions = base.boundary_condition,
        source_terms = source_terms,
        uEltype = eltype(p),
    )
    return semi, turbines
end

# function build_semi(p_list, base)
#     # turbines
#     # turbines = [Turbine(; x0 = p[1], y0 = p[2], u_rated, u_in, u_out, h_min) for p in p_list]

#     # source terms
#     # sponge_source = make_sponge_source(; Lx, σ_max=σₘ)
#     friction_source = source_terms_manning_friction
#     # turbine_source = make_turbine_source(turbines)
#     # all_source_terms = combine_source_terms([friction_source, turbine_source])
#     all_source_terms = combine_source_terms([friction_source])

#     semi = SemidiscretizationHyperbolic(
#         base.mesh, base.equations, base.initial_condition, base.solver;
#         boundary_conditions = base.boundary_condition,
#         source_terms = all_source_terms,
#     )
#     return semi
# end

function objective(p, base; tspan, saveat = 0.05, rho = 1000.0*20^3)
    semi, turbines = build_semi(p, base)
    ode = semidiscretize(semi, tspan)

    # Promote u0 to eltype(p) so ForwardDiff Dual numbers propagate through
    # the ODE state.  When p is plain Float64 this is a no-op.
    T = eltype(p)
    ode = remake(ode; u0 = T.(ode.u0))

    alive_callback = AliveCallback(alive_interval = 100)
    # summary_callback = SummaryCallback()
    callbacks = CallbackSet(alive_callback)

    stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

    sol = solve(ode, SSPRK43(stage_limiter!);
                dt = Δt,
                ode_default_options()...,
                callback = callbacks,
                adaptive = false,
                saveat)

    E, _ = compute_total_turbine_energy(sol, semi, turbines; rho = rho)
    # E = 0.0
    return -E
end

base = build_base_simulation()
p = [Lx / 6, 0.2 * Ly]  # flat vector: [x1, y1, ...] for ForwardDiff compatibility

# J0 = objective(p, base; tspan = (0.0, 0.25 * T), rho = 1000.0*D^3)

# println("J0 = objective(p) = ", J0)

# g = ForwardDiff.gradient(p -> objective(p, base; tspan = (0.0, T), rho = 1000.0*D^3), p)
# println("gradient at p0 = ", g)

# -------------------------------------------------------------------
# Gradient descent for optimal turbine placement
# -------------------------------------------------------------------
# Objective is -E (negative total energy), so minimizing it maximizes energy.
# Each evaluation traces ForwardDiff through the ODE solve, so iterations
# are expensive — keep max_iter modest and tune α / tspan_opt to taste.

function clip_to_domain!(p;
                         margin = 2 * 0.5,
                         headland_center = (Lx/2, Ly),
                         headland_radius = Rₛ)
    # Two constraints (both expanded by `margin` so the bump support stays clear):
    #  (1) Box: turbine inside [margin, Lx-margin] × [margin, Ly-margin].
    #  (2) Headland: turbine outside disk of radius (headland_radius + margin)
    #      centered on `headland_center` (the shallow shoal at the top).
    n = length(p) ÷ 2
    cx, cy = headland_center
    rmin   = headland_radius + margin

    for i in 1:n
        # 1) Box clip
        p[2i-1] = clamp(p[2i-1], margin, Lx - margin)
        p[2i]   = clamp(p[2i],   margin, Ly - margin)

        # 2) Push radially out of the headland disk if inside
        dx = p[2i-1] - cx
        dy = p[2i]   - cy
        d2 = dx^2 + dy^2
        if d2 < rmin^2
            d = sqrt(d2)
            if d < 1e-12
                # Degenerate: turbine right at headland center — pick a default direction.
                p[2i-1] = cx
                p[2i]   = cy - rmin
            else
                scale = rmin / d
                p[2i-1] = cx + dx * scale
                p[2i]   = cy + dy * scale
            end
            # Re-clip in case the radial projection pushed us outside the box
            p[2i-1] = clamp(p[2i-1], margin, Lx - margin)
            p[2i]   = clamp(p[2i],   margin, Ly - margin)
        end
    end
    return p
end

function gradient_descent(p0, base;
                          tspan_opt = (0.0, 0.25 * T),
                          rho = 1000.0 * D^3,
                          α = 0.5,             # step size in domain units
                          max_iter = 20,
                          grad_tol = 1e-3,
                          normalize_step = true,
                          csv_path = joinpath(output_directory, "gd_history.csv"),
                          run_label = "p0=$(p0)")
    p = copy(p0)
    f = p -> objective(p, base; tspan = tspan_opt, rho = rho)
    history = Tuple{Int, Float64, Vector{Float64}}[]   # (iter, J, p)

    n = length(p0)
    mkpath(dirname(csv_path))
    # First call (file missing): create + write header. Later calls: append only.
    if !isfile(csv_path)
        open(csv_path, "w") do io
            coords = join(["x$i,y$i" for i in 1:(n÷2)], ",")
            grads  = join(["gx$i,gy$i" for i in 1:(n÷2)], ",")
            println(io, "iter,J,E,grad_norm,$coords,$grads")
        end
    end
    # Title line separating this run from previous ones (prefixed with `#`
    # so CSV parsers with `comment="#"` skip it cleanly).
    open(csv_path, "a") do io
        println(io, "# === run: $(run_label)  tspan=$(tspan_opt)  α=$α ===")
    end

    for k in 0:max_iter
        result = ForwardDiff.DiffResults.GradientResult(p)
        ForwardDiff.gradient!(result, f, p)
        J = ForwardDiff.DiffResults.value(result)
        ∇J = ForwardDiff.DiffResults.gradient(result)
        gnorm = sqrt(sum(abs2, ∇J))

        push!(history, (k, J, copy(p)))
        @info "iter=$k  J=$(J)  E=$(-J)  ‖∇J‖=$(gnorm)  p=$(p)"

        # Append one CSV row, flushed each iteration so partial runs are recoverable
        open(csv_path, "a") do io
            row = string(k, ",", J, ",", -J, ",", gnorm, ",",
                         join(p, ","), ",", join(∇J, ","))
            println(io, row)
        end

        if gnorm < grad_tol
            @info "Gradient below tolerance — stopping"
            break
        end
        k == max_iter && break

        # Normalized step keeps each iteration moving by ~α regardless of |∇J|.
        # step = normalize_step ? α * ∇J / gnorm : α * ∇J
        # p .-= step
        # clip_to_domain!(p)

        d = ∇J ./ gnorm
        αk = α
        p_try = similar(p)                                                                                                                                  
        accepted = false
        while αk >= 1e-5                                                                                                                                    
            p_try .= p .- αk .* d                                                                                                                           
            clip_to_domain!(p_try)
            J_try = f(p_try)                                                                                                                                
            if J_try < J                                                                                                                                    
                p .= p_try
                accepted = true                                                                                                                             
                break                                             
            end
            αk *= 0.5
        end
        accepted || (@info "step rejected, stopping"; break)

    end
    return p, history
end

# Multi-start: try several initial points; everything goes into one CSV
p0_list = [
    # [Lx/6,    0.2 * Ly], 
    [Lx/2,    0.25 * Ly]
    # [Lx/2 - 6, Ly/4],     
    # [Lx/2 + 6, Ly/4],   
    # [Lx/2,    Ly/3],       
    # [3Lx/4,   Ly/4],       
]

# best_p = zero(p0_list)
# best_E = 0.0
for p0 in p0_list
    p_star, hist = gradient_descent(copy(p0), base;
                                    tspan_opt = (0.0, 0.25 * T),
                                    α = 0.5,
                                    max_iter = 200,
                                    grad_tol = 1e-5,
                                    run_label = "from p0=$(p0)")
    E_star = -hist[end][2]
    @info "from $p0  →  $p_star  E = $E_star"
    # if E_star > best_E
    #     best_E = E_star
    #     best_p = p_star
    # end
end

# println("\nBEST: optimal p = ", best_p, "   E = ", best_E)

# -------------------------------------------------------------------
# L-BFGS via Optim.jl (recommended over plain GD)
# -------------------------------------------------------------------
# Use Optim.only_fg!: when Optim wants both value and gradient, do a
# single ForwardDiff sweep; when it only wants a value (line search),
# do a plain Float64 forward solve — much cheaper than AD.

# using Optim

# function lbfgs_optimize(p0, base;
#                         tspan_opt = (0.0, 0.25 * T),
#                         rho = 1000.0 * D^3,
#                         margin = 2 * 0.5,
#                         iterations = 30,
#                         g_tol = 1e-6,
#                         show_trace = true)

#     obj    = p -> objective(p, base; tspan = tspan_opt, rho = rho)
#     n_eval = Ref(0)
#     n_grad = Ref(0)

#     function f_cb(p)
#         n_eval[] += 1
#         J = obj(p)                  # plain Float64 — used by line search
#         @info "[ feval #$(n_eval[])]  J=$J  p=$p"
#         return J
#     end

#     function g!_cb(G, p)
#         n_grad[] += 1
#         ForwardDiff.gradient!(G, obj, p)
#         @info "[grad  #$(n_grad[])]  ‖∇J‖=$(sqrt(sum(abs2, G)))  p=$p"
#         return G
#     end

#     # Box constraints: keep bump support inside the domain
#     n_turb = length(p0) ÷ 2
#     lo = repeat([margin,    margin   ], n_turb)
#     hi = repeat([Lx-margin, Ly-margin], n_turb)

#     res = optimize(f_cb, g!_cb, lo, hi, copy(p0),
#                    Fminbox(LBFGS()),
#                    Optim.Options(iterations  = iterations,
#                                  g_tol       = g_tol,
#                                  show_trace  = show_trace,
#                                  show_every  = 1))

#     return Optim.minimizer(res), Optim.minimum(res), res
# end

# p_star, J_star, res = lbfgs_optimize(p, base;
#                                      tspan_opt = (0.0, 0.25 * T),
#                                      iterations = 30)

# println("\noptimal p = ", p_star)
# println("final energy E = ", -J_star)
# println(res)