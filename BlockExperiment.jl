using Gmsh
using Gridap
using GridapGmsh
using GridapPETSc
using GridapMakie
using GLMakie

# PART 1: MESHING
gmsh.initialize()
gmsh.model.add("box2d")

# 2D rectangle: 10 wide, 20 tall (flow in y direction)
rect = gmsh.model.occ.addRectangle(0, 0, 0, 10, 20)

# Small fin protruding from x=0 wall
fin = gmsh.model.occ.addRectangle(0, 8, 0, 2, 4)

out, _ = gmsh.model.occ.fragment([(2, rect)], [(2, fin)])
gmsh.model.occ.synchronize()

tol = 0.01
fluid_tags = []
fin_tags = []
for s in gmsh.model.getEntities(2)
    tag = s[2]
    xmin,ymin,_,xmax,ymax,_ = gmsh.model.getBoundingBox(2, tag)
    cx = (xmin+xmax)/2
    if cx > 2.0
        push!(fluid_tags, tag)
    else
        push!(fin_tags, tag)
    end
end

fluid_tags = [3]
fin_tags = [2]

gmsh.model.addPhysicalGroup(2, fluid_tags, 1)
gmsh.model.setPhysicalName(2, 1, "fluid")
gmsh.model.addPhysicalGroup(2, fin_tags, 2)
gmsh.model.setPhysicalName(2, 2, "fin")

inlet_tags = []
outlet_tags = []
no_slip_tags = []
flux_tags = []
interface_tags = []

for c in gmsh.model.getEntities(1)
    tag = c[2]
    xmin,ymin,_,xmax,ymax,_ = gmsh.model.getBoundingBox(1, tag)
    parents = gmsh.model.getAdjacencies(1, tag)[1]

    if ymax < tol && ymin > -tol
        push!(inlet_tags, tag)
    elseif ymin > 20-tol
        push!(outlet_tags, tag)
    elseif 2 in parents && 3 in parents
        push!(interface_tags, tag)
    elseif xmax < tol && 2 in parents
        push!(flux_tags, tag)
    elseif 3 in parents
        push!(no_slip_tags, tag)
    end
end

gmsh.model.addPhysicalGroup(1, inlet_tags, 101)
gmsh.model.setPhysicalName(1, 101, "inlet")
gmsh.model.addPhysicalGroup(1, outlet_tags, 102)
gmsh.model.setPhysicalName(1, 102, "outlet")
gmsh.model.addPhysicalGroup(1, no_slip_tags, 103)
gmsh.model.setPhysicalName(1, 103, "no_slip")
gmsh.model.addPhysicalGroup(1, flux_tags, 104)
gmsh.model.setPhysicalName(1, 104, "flux")
gmsh.model.addPhysicalGroup(1, interface_tags, 105)
gmsh.model.setPhysicalName(1, 105, "interface")

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.5)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1.0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.mesh.generate(2)
gmsh.write("box2d.msh")
gmsh.finalize()

# PART 2: FE SPACES - FLOW
u_inlet = VectorValue(0.0, 1.0)
no_slip_vec = VectorValue(0.0, 0.0)

model = GmshDiscreteModel("box2d.msh")
Ω_f = Interior(model, tags="fluid")
Ω_s = Interior(model, tags="fin")
Γ = InterfaceTriangulation(Ω_s, Ω_f)

order = 2
reff_u = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
reff_p = ReferenceFE(lagrangian, Float64, order-1)

V = TestFESpace(Ω_f, reff_u, dirichlet_tags=["inlet","no_slip"])
U = TrialFESpace(V, [u_inlet, no_slip_vec])
Q = TestFESpace(Ω_f, reff_p, dirichlet_tags=["outlet"])
P = TrialFESpace(Q, [0.0])

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

println("Velocity DOFs: ", num_free_dofs(U))
println("Pressure DOFs: ", num_free_dofs(P))

# PART 3: WEAK FORM - FLOW
degree = 2*order
dΩ_f = Measure(Ω_f, degree)

μ = 1e-3

a((u,p),(v,q)) = ∫( μ*(∇(u)⊙∇(v)) - (∇⋅v)*p + q*(∇⋅u) )dΩ_f
l((v,q)) = ∫( 0.0*q )dΩ_f 

# PART 4: SOLVE STOKES EQUATION
options = [
    "-ksp_type", "preonly",
    "-pc_type", "lu",
    "-pc_factor_mat_solver_type", "mumps"
]

GridapPETSc.with(args=options) do
    println("Assembling Stokes...")
    op_stokes = AffineFEOperator(a, l, X, Y)
    println("Solving Stokes...")
    solver = PETScLinearSolver()
    uh, ph = solve(solver, op_stokes)
    println("Done Stokes")

    fig = Figure(size=(1000, 600))
    ax1 = Axis(fig[1,1], title="Velocity magnitude", aspect=DataAspect())
    ax2 = Axis(fig[1,2], title="Pressure", aspect=DataAspect())
    p1 = plot!(ax1, Ω_f, norm ∘ uh)
    p2 = plot!(ax2, Ω_f, ph)
    Colorbar(fig[2,1], p1, vertical=false, label="Velocity magnitude")
    Colorbar(fig[2,2], p2, vertical=false, label="Pressure")
    save("stokes_result.png", fig)
    println("Plot saved to stokes_result.png")



    #PART 5: SOLVE HEAT EQUATION
    # Material properties
    k_f = 0.6       # W/m·K water
    k_s = 200.0     # W/m·K aluminum fin
    ρ = 1000.0      # kg/m³
    cp = 4186.0     # J/kg·K
    T_inlet = 40.0 + 273.15
    q_flux = 1100.0  # W/m²

    # Measures
    dΩ_f = Measure(Ω_f, degree)
    dΩ_s = Measure(Ω_s, degree)
    dΓ = Measure(Γ, degree)
    dΓ_flux = Measure(BoundaryTriangulation(model, tags="flux"), degree)

    # Normal at interface (pointing from fluid to solid)
    n_Γ = get_normal_vector(Γ)

    # Nitsche penalty
    y = 1000.0

    # FE spaces for temperature
    reff_T = ReferenceFE(lagrangian, Float64, order)
    S_f = TestFESpace(Ω_f, reff_T, conformity=:H1, dirichlet_tags=["inlet"])
    S_s = TestFESpace(Ω_s, reff_T, conformity=:H1)
    Tf_trial = TrialFESpace(S_f, T_inlet)
    Ts_trial = TrialFESpace(S_s)

    Y_T = MultiFieldFESpace([S_f, S_s])
    X_T = MultiFieldFESpace([Tf_trial, Ts_trial])

    # Weak form
    a_T((tf,ts),(sf,ss)) =
        # Fluid diffusion + convection
        ∫( k_f*(∇(tf)⋅∇(sf)) + ρ*cp*((uh⋅∇(tf))*sf) )dΩ_f +
        # Solid conduction
        ∫( k_s*(∇(ts)⋅∇(ss)) )dΩ_s +
        # Nitsche interface coupling (solid=plus, fluid=minus)
        ∫( - k_s*(∇(ts).plus⋅n_Γ.plus)*ss.plus
        - k_s*(∇(ss).plus⋅n_Γ.plus)*ts.plus
        + k_f*(∇(tf).minus⋅n_Γ.minus)*sf.minus
        + k_f*(∇(sf).minus⋅n_Γ.minus)*tf.minus
        + k_s*(∇(ss).plus⋅n_Γ.plus)*tf.minus
        + k_f*(∇(sf).minus⋅n_Γ.minus)*ts.plus
        - k_s*(∇(ts).plus⋅n_Γ.plus)*sf.minus
        - k_f*(∇(tf).minus⋅n_Γ.minus)*ss.plus
        + y*(ts.plus - tf.minus)*(ss.plus - sf.minus) )dΓ

    l_T((sf,ss)) =
        ∫( q_flux * ss )dΓ_flux

    # Solve
    println("Assembling temperature...")
    op_T = AffineFEOperator(a_T, l_T, X_T, Y_T)
    Th_f, Th_s = solve(solver, op_T)
    println("Temperature done")
    

    # Plot
    fig2 = Figure(size=(1200, 500))
    ax1 = Axis(fig2[1,1], title="Fluid temperature", aspect=DataAspect())
    ax2 = Axis(fig2[1,2], title="Fin temperature", aspect=DataAspect())
    p1 = plot!(ax1, Ω_f, Th_f)
    p2 = plot!(ax2, Ω_s, Th_s)
    Colorbar(fig2[2,1], p1, vertical=false, label="T fluid (K)")
    Colorbar(fig2[2,2], p2, vertical=false, label="T solid (K)")
    save("temperature_result.png", fig2)
    println("Temperature plot saved")
end