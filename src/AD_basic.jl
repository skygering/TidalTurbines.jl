using HOHQMesh
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater
using TrixiBottomTopography
using CairoMakie
using Trixi2Vtk
using ForwardDiff

using LinearAlgebra
using Plots
using Test
##
equations = CompressibleEulerEquations2D(1.4)

solver = DGSEM(3, flux_central)
mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level = 2, n_cells_max = 10^5,
                periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_density_wave, solver;
                                    boundary_conditions = boundary_condition_periodic)

J = jacobian_ad_forward(semi);
size(J)
@test size(J) == (1024, 1024) #src

# Next, we compute the eigenvalues of the Jacobian.

λ = eigvals(J)
Plots.scatter(real.(λ), imag.(λ), label = "central flux")

# As you can see here, the maximal real part is close to zero.

relative_maximum = maximum(real, λ) / maximum(abs, λ)
@test 3.0e-10 < relative_maximum < 9.0e-10 #src

# Interestingly, if we add dissipation by switching to the `flux_lax_friedrichs`
# at the interfaces, the maximal real part of the eigenvalues increases.

solver = DGSEM(3, flux_lax_friedrichs)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_density_wave, solver;
                                    boundary_conditions = boundary_condition_periodic)

J = jacobian_ad_forward(semi)
λ = eigvals(J)

Plots.scatter!(real.(λ), imag.(λ), label = "Lax-Friedrichs flux")

##
equations = CompressibleEulerEquations1D(1.4)

solver = DGSEM(3, flux_central)
mesh = TreeMesh((-1.0,), (1.0,), initial_refinement_level = 6, n_cells_max = 10^5,
                periodicity = true)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_density_wave, solver;
                                    boundary_conditions = boundary_condition_periodic)

J = jacobian_ad_forward(semi)

λ = eigvals(J)

Plots.scatter(real.(λ), imag.(λ), label = "central flux")

## Differentiating through a complete simulation
Pkg.add("OrdinaryDiffEqLowOrderRK")
using OrdinaryDiffEqLowOrderRK

function energy_at_final_time(k) # k is the wave number of the initial condition
    equations = LinearScalarAdvectionEquation2D(1.0, -0.3)
    mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level = 3,
                    n_cells_max = 10^4, periodicity = true)
    solver = DGSEM(3, flux_lax_friedrichs)
    initial_condition = (x, t, equation) -> begin
        x_trans = Trixi.x_trans_periodic_2d(x - equation.advection_velocity * t)
        return SVector(sinpi(k * sum(x_trans)))
    end
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        uEltype = typeof(k),
                                        boundary_conditions = boundary_condition_periodic)
    ode = semidiscretize(semi, (0.0, 1.0))
    sol = solve(ode, FRK65(), dt = 0.05, adaptive = false, save_everystep = false)
    return Trixi.integrate(energy_total, sol.u[end], semi)
end

k_values = range(0.9, 1.1, length = 101)

Plots.plot(k_values, energy_at_final_time.(k_values), label = "Energy")

# You see a plot of a curve that resembles a parabola with local maximum around `k = 1.0`.
# Why's that? Well, the domain is fixed but the wave number changes. Thus, if the wave number is
# not chosen as an integer, the initial condition will not be a smooth periodic function in the
# given domain. Hence, the dissipative surface flux (`flux_lax_friedrichs` in this example)
# will introduce more dissipation. In particular, it will introduce more dissipation for "less smooth"
# initial data, corresponding to wave numbers `k` further away from integers.

# We can compute the discrete derivative of the energy at the final time with respect to the wave
# number `k` as follows.

round(ForwardDiff.derivative(energy_at_final_time, 1.0), sigdigits = 2)
@test round(ForwardDiff.derivative(energy_at_final_time, 1.0), sigdigits = 2) == 1.4e-5 #src

# This is rather small and we can treat it as zero in comparison to the value of this derivative at
# other wave numbers `k`.

dk_values = ForwardDiff.derivative.((energy_at_final_time,), k_values);

Plots.plot(k_values, dk_values, label = "Derivative")

## If you remember basic calculus, a sufficient condition for a local maximum is that the first derivative
# vanishes and the second derivative is negative. We can also check this discretely.

second_derivative = round(ForwardDiff.derivative(k -> Trixi.ForwardDiff.derivative(energy_at_final_time,
                                                                                   k), 1.0),
                          sigdigits = 2)
@test second_derivative ≈ -0.9 #src

# Having seen this application, let's break down what happens step by step.

function energy_at_final_time(k) # k is the wave number of the initial condition
    equations = LinearScalarAdvectionEquation2D(1.0, -0.3)
    mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level = 3,
                    n_cells_max = 10^4, periodicity = true)
    solver = DGSEM(3, flux_lax_friedrichs)
    initial_condition = (x, t, equation) -> begin
        x_trans = Trixi.x_trans_periodic_2d(x - equation.advection_velocity * t)
        return SVector(sinpi(k * sum(x_trans)))
    end
    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        uEltype = typeof(k),
                                        boundary_conditions = boundary_condition_periodic)
    ode = semidiscretize(semi, (0.0, 1.0))
    sol = solve(ode, FRK65(), dt = 0.05, adaptive = false, save_everystep = false)
    return Trixi.integrate(energy_total, sol.u[end], semi)
end

k = 1.0
round(ForwardDiff.derivative(energy_at_final_time, k), sigdigits = 2)
@test round(ForwardDiff.derivative(energy_at_final_time, 1.0), sigdigits = 2) == 1.4e-5 #src

# When calling `ForwardDiff.derivative(energy_at_final_time, k)` with `k=1.0`, ForwardDiff.jl
# will basically use the chain rule and known derivatives of existing basic functions
# to calculate the derivative of the energy at the final time with respect to the
# wave number `k` at `k0 = 1.0`. To do this, ForwardDiff.jl uses dual numbers, which
# basically store the result and its derivative w.r.t. a specified parameter at the
# same time. Thus, we need to make sure that we can treat these `ForwardDiff.Dual`
# numbers everywhere during the computation. Fortunately, generic Julia code usually
# supports these operations. The most basic problem for a developer is to ensure
# that all types are generic enough, in particular the ones of internal caches.

# The first step in this example creates some basic ingredients of our simulation.
equations = LinearScalarAdvectionEquation2D(1.0, -0.3)
mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level = 3, n_cells_max = 10^4,
                periodicity = true)
solver = DGSEM(3, flux_lax_friedrichs);

# These do not have internal caches storing intermediate values of the numerical
# solution, so we do not need to adapt them. In fact, we could also define them
# outside of `energy_at_final_time` (but would need to take care of globals or
# wrap everything in another function).

# Next, we define the initial condition
initial_condition = (x, t, equation) -> begin
    x_trans = Trixi.x_trans_periodic_2d(x - equation.advection_velocity * t)
    return SVector(sinpi(k * sum(x_trans)))
end;
# as a closure capturing the wave number `k` passed to `energy_at_final_time`.
# If you call `energy_at_final_time(1.0)`, `k` will be a `Float64`. Thus, the
# return values of `initial_condition` will be `SVector`s of `Float64`s. When
# calculating the `ForwardDiff.derivative`, `k` will be a `ForwardDiff.Dual` number.
# Hence, the `initial_condition` will return `SVector`s of `ForwardDiff.Dual`
# numbers.

# The semidiscretization `semi` uses some internal caches to avoid repeated allocations
# and speed up the computations, e.g. for numerical fluxes at interfaces. Thus, we
# need to tell Trixi.jl to allow `ForwardDiff.Dual` numbers in these caches. That's what
# the keyword argument `uEltype=typeof(k)` in
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    uEltype = typeof(k),
                                    boundary_conditions = boundary_condition_periodic)

# does. This is basically the only part where you need to modify your standard Trixi.jl
# code to enable automatic differentiation. From there on, the remaining steps
ode = semidiscretize(semi, (0.0, 1.0))
sol = solve(ode, FRK65(), dt = 0.05, adaptive = false, save_everystep = false)
round(Trixi.integrate(energy_total, sol.u[end], semi), sigdigits = 5)
@test round(Trixi.integrate(energy_total, sol.u[end], semi), sigdigits = 5) == 0.25 #src