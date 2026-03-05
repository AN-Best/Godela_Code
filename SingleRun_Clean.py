from OFOAM_Functions import create_mesh, run_solver, view_results, get_fin_volume
import os

case_dir      = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_case"
stl_coldplate = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_inst_stl/coldplate_inst/Body_6.stl"
stl_fins      = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_inst_stl/coldplate_inst/Body_3.stl"

T_inlet    = 313.15        # K  (40°C)
flow_rate  = 1.65e-3 / 60  # m³/s
heat_flux  = 1100           # W

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
    end_time      = 50,
    write_interval= 5,
    show_plots=True
)


# STEP 3: VIEW RESULTS
results = view_results(
    case_dir   = case_dir,
    T_inlet    = T_inlet,
    heat_flux  = heat_flux,
    show_plots = True
)

print(f"\n Final Result ")
print(f"T_max = {results['T_max']:.2f} K  ({results['T_max']-273.15:.2f}°C)")
print(f"R_th  = {results['R_th']:.6f} K/W")
print(f"Delta_P = {results['delta_p']:.2f} Pa")

