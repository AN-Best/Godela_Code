using Gmsh
using Gridap
using GridapGmsh
using GridapPETSc
using GridapMakie
using GLMakie

#PART 1: MESHING
plotter = true

if !isfile("coldplate.msh")

    gmsh.initialize()
    gmsh.model.add("coldplate")

    path = "/home/aaronbest/Desktop/GodelaChallenge/coldplate/coldplate_inst10.stp"
    v = gmsh.model.occ.import_shapes(path)
    gmsh.model.occ.synchronize()

    name_array = ["case","bottom_cover","fin","turret1","turret2","fluid"]
    for p in v
        dim = Int(p[1])
        tag = Int(p[2])
        gmsh.model.addPhysicalGroup(dim, [tag], tag)
        gmsh.model.setPhysicalName(dim, tag, name_array[tag])
    end

    flux_surface_tag = [258]
    inlet_surface_tag = [1499]
    outlet_surface_tag = [1508]
    fin_surface_tag = []
    no_slip_tag = []

    for surface in gmsh.model.getEntities(2)
        tag = surface[2]
        parent_volume = gmsh.model.getAdjacencies(2, tag)[1]
        if 3 in parent_volume && !(tag in flux_surface_tag)
            push!(fin_surface_tag, tag)
        end
        if 6 in parent_volume &&
           !(tag in inlet_surface_tag) &&
           !(tag in outlet_surface_tag) &&
           !(4 in parent_volume) &&
           !(5 in parent_volume)
            push!(no_slip_tag, tag)
        end
    end

    gmsh.model.addPhysicalGroup(2, flux_surface_tag, 100)
    gmsh.model.setPhysicalName(2, 100, "flux")
    gmsh.model.addPhysicalGroup(2, inlet_surface_tag, 101)
    gmsh.model.setPhysicalName(2, 101, "inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surface_tag, 102)
    gmsh.model.setPhysicalName(2, 102, "outlet")
    gmsh.model.addPhysicalGroup(2, fin_surface_tag, 103)
    gmsh.model.setPhysicalName(2, 103, "fin")
    gmsh.model.addPhysicalGroup(2, no_slip_tag, 104)
    gmsh.model.setPhysicalName(2, 104, "no_slip")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.01)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1.0)

    gmsh.model.mesh.generate(3)
    gmsh.write("coldplate.msh")

    if plotter; gmsh.fltk.run(); end
    gmsh.finalize()
end

#PART 2: FLOW
flowrate = (1.65e-3/60)*1e9
μ = 1e-3

model = GmshDiscreteModel("coldplate.msh")

order = 2
degree = 2*order

Ω_fluid = Interior(model, tags="fluid")
dΩ_fluid = Measure(Ω_fluid, degree)

reff_u = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)
reff_p = ReferenceFE(lagrangian, Float64, order-1)

# Geometric inlet: faces flat at z=77 (inlet end cap only)
Γ_all_bound_pre = BoundaryTriangulation(Ω_fluid)
mask_inlet_geom = map(get_cell_coordinates(Γ_all_bound_pre)) do cell_coords
    zs = [c[3] for c in cell_coords]
    minimum(zs) > 76.9 && maximum(zs) < 77.1
end
Γ_inlet_geom = Triangulation(Γ_all_bound_pre, findall(mask_inlet_geom))
dΓ_inlet_geom = Measure(Γ_inlet_geom, degree)
inlet_area_geom = sum(∫(1.0)dΓ_inlet_geom)
u_inlet_actual = flowrate / inlet_area_geom
println("Geometric inlet face count: ", num_cells(Γ_inlet_geom))
println("Geometric inlet area:       ", inlet_area_geom, " mm²")
println("Theoretical inlet area:     ", π*1.9^2, " mm²")
println("u_inlet_actual:             ", u_inlet_actual, " mm/s")

# Geometric outlet: faces flat at z=-12 (outlet end cap only)
mask_outlet_geom = map(get_cell_coordinates(Γ_all_bound_pre)) do cell_coords
    zs = [c[3] for c in cell_coords]
    minimum(zs) > -12.1 && maximum(zs) < -11.9
end
Γ_outlet_geom = Triangulation(Γ_all_bound_pre, findall(mask_outlet_geom))
dΓ_outlet_geom = Measure(Γ_outlet_geom, degree)
outlet_area_geom = sum(∫(1.0)dΓ_outlet_geom)
println("Geometric outlet face count: ", num_cells(Γ_outlet_geom))
println("Geometric outlet area:       ", outlet_area_geom, " mm²")

# Penalty method to weakly enforce inlet velocity on geometric faces
γ_inlet = 1e6 * μ
u_in = VectorValue(0.0, 0.0, -u_inlet_actual)

V = TestFESpace(Ω_fluid, reff_u, conformity=:H1, dirichlet_tags=["no_slip"])
U = TrialFESpace(V, VectorValue(0.0,0.0,0.0))
Q = TestFESpace(Ω_fluid, reff_p, conformity=:H1, dirichlet_tags=["outlet"])
P = TrialFESpace(Q, 0.0)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

a((u,p),(v,q)) = ∫( μ*(∇(u)⊙∇(v)) - (∇⋅v)*p + q*(∇⋅u) )dΩ_fluid +
                 ∫( γ_inlet*(u⋅v) )dΓ_inlet_geom

l((v,q)) = ∫( 0.0*q )dΩ_fluid +
           ∫( γ_inlet*(u_in⋅v) )dΓ_inlet_geom

options = [
    "-ksp_type", "preonly",
    "-pc_type", "lu",
    "-pc_factor_mat_solver_type", "mumps"
]

println("Velocity DOFs: ", num_free_dofs(U))
println("Pressure DOFs: ", num_free_dofs(P))

uh_global = nothing
ph_global = nothing

GridapPETSc.with(args=options) do
    println("Assembling Stokes...")
    op = AffineFEOperator(a, l, X, Y)
    solver = PETScLinearSolver()
    uh, ph = solve(solver, op)
    println("Max velocity: ", maximum(abs.(get_free_dof_values(uh))))
    println("Max pressure: ", maximum(abs.(get_free_dof_values(ph))))

    dΓi = Measure(Γ_inlet_geom, degree)
    dΓo = Measure(Γ_outlet_geom, degree)
    Γ_all_f = BoundaryTriangulation(Ω_fluid)
    dΓa = Measure(Γ_all_f, degree)
    println("Inlet flux:  ", sum(∫(uh⋅get_normal_vector(Γ_inlet_geom))dΓi), " mm³/s")
    println("Outlet flux: ", sum(∫(uh⋅get_normal_vector(Γ_outlet_geom))dΓo), " mm³/s")
    println("Total flux:  ", sum(∫(uh⋅get_normal_vector(Γ_all_f))dΓa), " mm³/s")
    println("Expected:    ", flowrate, " mm³/s")

    writevtk(Ω_fluid, "stokes_results", cellfields=["u"=>uh, "p"=>ph])
    global uh_global = uh
    global ph_global = ph

    fig = Figure(size=(1200, 500))
    ax1 = Axis3(fig[1,1], title="Velocity magnitude")
    ax2 = Axis3(fig[1,2], title="Pressure")
    p1 = plot!(ax1, Ω_fluid, norm ∘ uh)
    p2 = plot!(ax2, Ω_fluid, ph)
    Colorbar(fig[2,1], p1, vertical=false, label="Velocity (mm/s)")
    Colorbar(fig[2,2], p2, vertical=false, label="Pressure (Pa)")
    save("coldplate_stokes.png", fig)
    println("Stokes plot saved")
end

#PART 3 - FLUID TEMPERATURE WITH VOLUMETRIC HEAT SOURCE

k_f = 0.6e-3
ρ = 1.0e-6
cp = 4186.0
T_inlet = 40.0 + 273.15
geom_tol = 0.5

dΩ_f = Measure(Ω_fluid, degree)

mask_fluid_fin = map(get_cell_coordinates(Ω_fluid)) do cell_coords
    ys = [c[2] for c in cell_coords]
    maximum(ys) < 1.19 + geom_tol && minimum(ys) > -3.31 - geom_tol
end
Ω_fluid_fin = Triangulation(Ω_fluid, findall(mask_fluid_fin))
dΩ_fluid_fin = Measure(Ω_fluid_fin, degree)

fin_fluid_volume = sum(∫(1.0)dΩ_fluid_fin)
Q_vol = 1100.0 / fin_fluid_volume
println("Fluid volume near fins: ", fin_fluid_volume, " mm³")
println("Q_vol = ", Q_vol, " W/mm³")
println("Total heat input check: ", sum(∫(Q_vol * 1.0)dΩ_fluid_fin), " W")

reff_T = ReferenceFE(lagrangian, Float64, order)
S_f = TestFESpace(Ω_fluid, reff_T, conformity=:H1, dirichlet_tags=["inlet"])
Tf_trial = TrialFESpace(S_f, T_inlet)

h_mesh = 0.5

τ = Operation(x -> begin
    u_norm = max(norm(x), 1e-10)
    pe_local = ρ * cp * u_norm * h_mesh / (2 * k_f)
    if pe_local > 1
        h_mesh / (2 * ρ * cp * u_norm) * (1 - 1/pe_local)
    else
        h_mesh^2 / (6 * k_f)
    end
end)(uh_global)

a_f(tf, sf) =
    ∫( k_f*(∇(tf)⋅∇(sf)) + ρ*cp*((uh_global⋅∇(tf))*sf) + τ*(uh_global⋅∇(sf))*(uh_global⋅∇(tf)) )dΩ_f

l_f(sf) =
    ∫( Q_vol * sf + τ*(uh_global⋅∇(sf))*Q_vol )dΩ_fluid_fin

println("Temperature DOFs: ", num_free_dofs(Tf_trial))

GridapPETSc.with(args=options) do
    println("Assembling temperature...")
    op_f = AffineFEOperator(a_f, l_f, Tf_trial, S_f)
    solver = PETScLinearSolver()
    Th_f = solve(solver, op_f)
    println("Temperature done")
    println("Max fluid temp: ", maximum(get_free_dof_values(Th_f)), " K")
    println("Min fluid temp: ", minimum(get_free_dof_values(Th_f)), " K")

    dΓo_t = Measure(Γ_outlet_geom, degree)
    dΓi_t = Measure(Γ_inlet_geom, degree)
    heat_out = sum(∫( ρ*cp*(uh_global⋅get_normal_vector(Γ_outlet_geom))*Th_f )dΓo_t)
    heat_in  = sum(∫( ρ*cp*(uh_global⋅get_normal_vector(Γ_inlet_geom))*Th_f )dΓi_t)
    println("Net convective heat transport: ", heat_out + heat_in, " W")
    println("Expected: 1100 W")

    writevtk(Ω_fluid, "temperature_fluid", cellfields=["T"=>Th_f])

    fig = Figure(size=(800, 500))
    ax = Axis3(fig[1,1], title="Fluid temperature")
    p = plot!(ax, Ω_fluid, Th_f)
    Colorbar(fig[2,1], p, vertical=false, label="T (K)")
    save("coldplate_temperature.png", fig)
    println("Temperature plot saved")
end