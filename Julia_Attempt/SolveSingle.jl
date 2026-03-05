using Gridap
using GridapGmsh
using LinearAlgebra
using IterativeSolvers
using IncompleteLU

plotter = true

# ── Load mesh ─────────────────────────────────────────────────────────────────
model = GmshDiscreteModel("Julia_Attempt/coldplate.msh")

Ω_all     = Triangulation(model)
reffe_chk = ReferenceFE(lagrangian, Float64, 1)
V_chk     = TestFESpace(model, reffe_chk)

println("\n── Mesh Statistics ─────────────────────")
println("Total cells:  ", num_cells(Ω_all))
println("Total nodes:  ", num_free_dofs(V_chk))
println("────────────────────────────────────────")

# ── Physical parameters (mm units) ────────────────────────────────────────────
flowrate = 1.65e-3/60 * 1e9   # mm³/s
r        = 3.8/2               # mm
A_inlet  = π * r^2             # mm²
u_mag    = flowrate / A_inlet  # mm/s

μ        = 1.002e-3 * 1e-3    # kg/(mm·s)
ρ        = 998.0   * 1e-9     # kg/mm³
Cp       = 4182.0  * 1e6      # mJ/(kg·K)
κ_f      = 0.598   * 1e-3     # mW/(mm·K)
T_inlet  = 313.15              # K
power    = 1100.0  * 1e3      # mW

println("Inlet velocity: ", round(u_mag, digits=4), " mm/s")

# ── Porous medium parameters ───────────────────────────────────────────────────
fin_gap   = 0.23
fin_t     = 0.2
H_fin     = 3.0
fin_pitch = fin_t + fin_gap

κ_solid   = 400.0 * 1e-3
ε         = fin_gap / fin_pitch
κ_eff     = (fin_t*κ_solid + fin_gap*κ_f) / fin_pitch
D_h       = 2 * fin_gap * H_fin / (fin_gap + H_fin)
f_Re      = 96.0
K         = D_h^2 / (2 * f_Re)

println("Porosity ε:     ", round(ε,    digits=4))
println("Permeability K: ", round(K,    digits=6), " mm²")
println("κ_eff:          ", round(κ_eff,digits=4), " mW/(mm·K)")
println("μ/K:            ", round(μ/K,  digits=4))

# ── Triangulations ────────────────────────────────────────────────────────────
degree   = 1
Ω_free   = Triangulation(model, tags=["fluid"])
Ω_porous = Triangulation(model, tags=["fins"])
Ω_h      = Triangulation(model, tags=["fluid","fins"])

dΩ_free   = Measure(Ω_free,   degree)
dΩ_porous = Measure(Ω_porous, degree)
dΩ        = Measure(Ω_h,      degree)

# ── FE spaces: P1/P1 + pressure stabilization (PSPG) ─────────────────────────
# P1/P1 is used instead of Taylor-Hood to keep DOF count manageable on 16GB RAM.
# PSPG stabilization restores inf-sup stability.
reffe_u = ReferenceFE(lagrangian, VectorValue{3,Float64}, 1)
reffe_p = ReferenceFE(lagrangian, Float64, 1)

V = TestFESpace(model, reffe_u, conformity=:H1,
                dirichlet_tags=["inlet","walls"],
                dirichlet_masks=[(true,true,true),(true,true,true)])
Q = TestFESpace(model, reffe_p, conformity=:H1,
                dirichlet_tags=["outlet"])
Y = MultiFieldFESpace([V, Q])

u_inlet = VectorValue(0.0, 0.0, -u_mag)
u_wall  = VectorValue(0.0, 0.0,  0.0)
U = TrialFESpace(V, [u_inlet, u_wall])
P = TrialFESpace(Q, 0.0)
X = MultiFieldFESpace([U, P])

println("\nFree DOFs — U: ", num_free_dofs(U), "  P: ", num_free_dofs(P))

# PSPG stabilization parameter τ_pspg = h²/(4μ) element-wise
# Here we use a global estimate with h ~ D_h
D_h    = 2 * 0.23 * 3.0 / (0.23 + 3.0)
τ_pspg = D_h^2 / (4 * μ)

f = VectorValue(0.0, 0.0, 0.0)

# Stokes + Darcy weak form with PSPG pressure stabilization
a_flow((u,p),(v,q)) =
    ∫( μ*(∇(v)⊙∇(u)) - (∇⋅v)*p + q*(∇⋅u)
       + τ_pspg*(∇(q)⋅∇(p)) )dΩ_free   +
    ∫( (μ/K)*(v⋅u)   - (∇⋅v)*p + q*(∇⋅u)
       + τ_pspg*(∇(q)⋅∇(p)) )dΩ_porous

l_flow((v,q)) =
    ∫( v⋅f )dΩ_free +
    ∫( v⋅f )dΩ_porous

op_flow = AffineFEOperator(a_flow, l_flow, X, Y)

A = op_flow.op.matrix
b = op_flow.op.vector

println("\nAssembling preconditioner...")
P_ilu = ilu(A, τ = 0.01)

println("Solving flow with GMRES...")
x, history = gmres(A, b; Pl=P_ilu, restart=50, maxiter=500, reltol=1e-8, log=true)
println("  Converged: ", history.isconverged, "  iters: ", history.iters)

xh = FEFunction(X, x)
uh, ph = xh

println("\nFlow solved!")
println("Max |u|: ", maximum(norm.(uh.free_values)), " mm/s")
println("Max |p|: ", maximum(abs.(ph.free_values)),  " MPa")

# ── Interface heat flux ────────────────────────────────────────────────────────
Γ           = BoundaryTriangulation(model, tags=["fluid_porous_interface"])
dΓ          = Measure(Γ, degree)
A_interface = sum(∫(1.0)dΓ)
q_flux      = power / A_interface
println("\nInterface area: ", round(A_interface, digits=2), " mm²")
println("Heat flux:      ", round(q_flux,       digits=4), " mW/mm²")

# ── Temperature spaces ────────────────────────────────────────────────────────
reffe_T = ReferenceFE(lagrangian, Float64, 1)
W       = TestFESpace(model, reffe_T, conformity=:H1, dirichlet_tags=["inlet"])
T_trial = TrialFESpace(W, T_inlet)

h_elem = D_h
u_ref  = max(u_mag, 1e-10)
τ_supg = h_elem / (2 * u_ref)

a_T(t, w) =
    ∫( ρ*Cp*(uh⋅∇(t))*w + κ_f  *(∇(t)⋅∇(w))
       + τ_supg*ρ*Cp*(uh⋅∇(w))*(uh⋅∇(t)) )dΩ_free   +
    ∫( ρ*Cp*(uh⋅∇(t))*w + κ_eff*(∇(t)⋅∇(w))
       + τ_supg*ρ*Cp*(uh⋅∇(w))*(uh⋅∇(t)) )dΩ_porous

l_T(w) = ∫( q_flux * w )dΓ

op_T = AffineFEOperator(a_T, l_T, T_trial, W)

A_T = op_T.op.matrix
b_T = op_T.op.vector

println("\nSolving temperature with GMRES...")
P_T   = ilu(A_T, τ = 0.01)
x_T, hist_T = gmres(A_T, b_T; Pl=P_T, restart=50, maxiter=500, reltol=1e-8, log=true)
println("  Converged: ", hist_T.isconverged, "  iters: ", hist_T.iters)

Th = FEFunction(T_trial, x_T)

println("\nTemperature solved!")
println("T_max: ", round(maximum(Th.free_values) - 273.15, digits=2), " °C")
println("T_min: ", round(minimum(Th.free_values) - 273.15, digits=2), " °C")
println("ΔT:    ", round(maximum(Th.free_values) - T_inlet, digits=4), " K")

# ── Key results ───────────────────────────────────────────────────────────────
ΔT   = maximum(Th.free_values) - T_inlet
R_th = ΔT / (power * 1e-3)
println("\n── Key Results ──────────────────────────")
println("Thermal resistance R_th: ", round(R_th, digits=4), " K/W")
println("────────────────────────────────────────")

writevtk(Ω_h, "fluid_results", cellfields=["u"=>uh, "p"=>ph, "T"=>Th])
println("\nResults saved to fluid_results.vtu")