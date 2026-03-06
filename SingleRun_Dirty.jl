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


    #Separate volumes into physical groups
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


    surfaces = gmsh.model.getEntities(2)  # Get all 2D entities (surfaces)
    println("Available surfaces:")
    for surface in surfaces
        tag = surface[2]
        bbox = gmsh.model.getBoundingBox(2, tag)
        parent_volume = gmsh.model.getAdjacencies(2,tag)[1]
        if 3 in parent_volume && !(tag in flux_surface_tag)
            push!(fin_surface_tag,tag)
        end
        if 6 in parent_volume
            push!(no_slip_tag,tag)
        end

    end

    #Bottom of pedestal flux
    gmsh.model.addPhysicalGroup(2,flux_surface_tag,100)
    gmsh.model.setPhysicalName(2,100,"flux")
    #Inlet of flow
    gmsh.model.addPhysicalGroup(2,inlet_surface_tag,101)
    gmsh.model.setPhysicalName(2,101,"inlet")
    #Outlet
    gmsh.model.addPhysicalGroup(2,outlet_surface_tag,102)
    gmsh.model.setPhysicalName(2,102,"outlet")
    #Fins
    gmsh.model.addPhysicalGroup(2,fin_surface_tag,103)
    gmsh.model.setPhysicalName(2,103,"fin")
    #No Slip
    gmsh.model.addPhysicalGroup(2,no_slip_tag,104)
    gmsh.model.setPhysicalName(2,104,"no_slip")


    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3.0)

    gmsh.model.mesh.generate(3)
    gmsh.write("coldplate.msh")


    if plotter
        gmsh.fltk.run() 
    end

    gmsh.finalize()

end

#PART 2: FLOW

#Parameters
flowrate = (1.65e-3/60)*1e9
r = 3.8/2.0
A = π*r^2
u_mag = flowrate/A
μ = 1e-3  

println("u_mag = ", u_mag)

u_inlet = VectorValue(0.0,0.0,-u_mag)
no_slip_vec = VectorValue(0.0,0.0,0.0)

model = GmshDiscreteModel("coldplate.msh")

order = 2

degree = 2*order
Ω_fluid = Interior(model, tags="fluid") 
dΩ_fluid = Measure(Ω_fluid, degree)

reff_u = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
reff_p = ReferenceFE(lagrangian,Float64,order-1)

V = TestFESpace(Ω_fluid, reff_u, conformity=:H1, dirichlet_tags=["inlet","no_slip"])
U = TrialFESpace(V, [u_inlet, no_slip_vec])
Q = TestFESpace(Ω_fluid, reff_p, conformity=:H1, dirichlet_tags=["outlet"])
P = TrialFESpace(Q, [0.0])

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

#Weak Form
a((u,p),(v,q)) = ∫( μ*(∇(u)⊙∇(v)) - (∇⋅v)*p + q*(∇⋅u) )dΩ_fluid
l((v,q)) = ∫( 0.0*q )dΩ_fluid

# Assemble and solve
options = [
    "-ksp_type", "gmres",
    "-ksp_rtol", "1.0e-6",
    "-ksp_max_it", "1000",
    "-pc_type", "jacobi",
    "-ksp_monitor"
]

println("Velocity DOFs: ", num_free_dofs(U))
println("Pressure DOFs: ", num_free_dofs(P))

GridapPETSc.with(args=options) do


    println("Starting assembly...")
    op = AffineFEOperator(a, l, X, Y)
    println("Assembly done, starting solve...")
    solver = PETScLinearSolver()
    uh, ph = solve(solver, op)
    println("Solve done")
    println("Max velocity: ", maximum(abs.(get_free_dof_values(uh))))
    println("Max pressure: ", maximum(abs.(get_free_dof_values(ph))))
    writevtk(Ω_fluid, "stokes_results", cellfields=["u"=>uh, "p"=>ph])

    fig = Figure(size=(1200, 500))
    ax1 = Axis3(fig[1,1], title="Velocity magnitude")
    ax2 = Axis3(fig[1,2], title="Pressure")
    p1 = plot!(ax1, Ω_fluid, norm ∘ uh)
    p2 = plot!(ax2, Ω_fluid, ph)
    Colorbar(fig[2,1], p1, vertical=false, label="Velocity (mm/s)")
    Colorbar(fig[2,2], p2, vertical=false, label="Pressure (Pa)")
    save("coldplate_stokes.png", fig)
    println("Plot saved")
end
