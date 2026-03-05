from OFOAM_Functions import create_mesh, run_solver, view_results, get_fin_volume
import os

case_dir      = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_case"
stl_coldplate = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_inst_stl/coldplate_inst/Body_6.stl"
stl_fins      = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_inst_stl/coldplate_inst/Body_3.stl"

base = "/home/aaronbest/Desktop/GodelaChallenge"


T_inlet    = 313.15        
flow_rate  = 1.65e-3 / 60  
heat_flux  = 1100           

#Number of iterations to run
N_runs = 3
for i in range(N_runs):

    case_folder = "Case" + str(i)
    stl_folder = "coldplate_inst" + str(i)

    case_dir = os.path.join(base,case_folder)
    stl_coldplate = os.path.join(base,"coldplate_inst_stl",stl_folder,"Body_6.stl")
    stl_fins = os.path.join(base,"coldplate_inst_stl",stl_folder,"Body_3.stl")

    # STEP 1: MESH
    mesh_exists = os.path.exists(f"{case_dir}/constant/solid/polyMesh/boundary")
    if mesh_exists:
        print("Mesh already exists — skipping meshing step.")
    else:
        print("No mesh found — creating mesh...")
        create_mesh(
            case_dir      = case_dir,
            stl_coldplate = stl_coldplate,
            stl_fins      = stl_fins,
            pad           = 0.005,
            nx=50, ny=50, nz=50
        )


    # STEP 2: SOLVE
    fin_volume = get_fin_volume(case_dir=case_dir)

    run_solver(
        case_dir      = case_dir,
        T_inlet       = T_inlet,
        flow_rate     = flow_rate,
        heat_flux     = heat_flux,
        fin_volume    = fin_volume,
        end_time      = 2500,
        write_interval= 50,
        show_plots=False
    )


results = view_results(
    case_dir   = os.path.join(base,"Case2"),
    T_inlet    = T_inlet,
    heat_flux  = heat_flux,
    show_plots = True
)