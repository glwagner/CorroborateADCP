using Oceananigans

# Domain
Nz = Nh = 64
Lz = 32
Lh = 2Lz

grid = RectilinearGrid(size=(Nh, Nh, Nz), extent=(Lh, Lh, Lz), halo=(3, 3, 3))

# Surface momentum flux
t★ = 10minutes # time-scale for increasing momentum flux
u★ = 1e-2 # [m s⁻¹] friction velocity
@inline x_momentum_flux(x, y, t, p) = - p.u★^2 * min(1, sqrt(t / p.t★))
u_top_bc = FluxBoundaryCondition(x_momentum_flux, parameters=(; u★, t★))
u_boundary_conditions = FieldBoundaryConditions(top=u_top_bc)
boundary_conditions = (; u=u_boundary_conditions)

# Surface waves
La = 0.1
const k = 2π / 10 # [m⁻¹] effective peak wavenumber
const uˢ = u★ / La # surface Stokes drift
@inline ∂z_uˢ(z, t) = k * uˢ * exp(2k * z)
stokes_drift = UniformStokesDrift(; ∂z_uˢ)

coriolis = nothing #FPlane(f=1e-4)
model = NonhydrostaticModel(; grid, coriolis, boundary_conditions, stokes_drift, advection=WENO5())

ϵ(x, y, z) = 1e-3 * u★ * rand()
set!(model, u=ϵ, v=ϵ, w=ϵ)

simulation = Simulation(model, Δt=1.0, stop_time=1hour)

progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:p] = Callback(progress, IterationInterval(10))

run!(simulation)


