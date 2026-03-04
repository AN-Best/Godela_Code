import os
import pyvista as pv
import subprocess
import shutil
import glob
import math
import matplotlib.pyplot as plt
import re

## CONFIGURATION
OPENFOAM_ENV = os.environ.copy()
OPENFOAM_ENV["PATH"] = "/opt/openfoam11/bin:" + OPENFOAM_ENV["PATH"]
OPENFOAM_ENV["WM_PROJECT_DIR"] = "/opt/openfoam11"

def run_cmd(cmd, cwd):
    print(f"\nRunning: {cmd}")
    result = subprocess.run(
        cmd, shell=True, cwd=cwd,
        env=OPENFOAM_ENV, executable="/bin/bash"
    )
    if result.returncode != 0:
        print(f"ERROR: {cmd} failed with return code {result.returncode}")
        exit(1)
    print(f"Done: {cmd}")

def write_file(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    print(f"Written: {path}")

## Fin Volume

def get_fin_volume(case_dir):
    """
    Runs checkMesh on region1 and extracts the total cell volume.
    Must be called after meshing.
    """
    result = subprocess.run(
        "source /opt/openfoam11/etc/bashrc && checkMesh -region region1",
        shell=True, cwd=case_dir, env=OPENFOAM_ENV,
        executable="/bin/bash", capture_output=True, text=True
    )

    for line in result.stdout.splitlines():
        if "Total volume" in line:
            match = re.search(r"Total volume = ([\d.e+-]+)\.", line)
            if match:
                volume = float(match.group(1))
                print(f"Fin volume: {volume:.6e} m³")
                return volume

    raise RuntimeError("Could not extract fin volume from checkMesh output")

## CREATE MESH

def create_mesh(case_dir, stl_coldplate, stl_fins, pad=0.005, nx=50, ny=50, nz=50, force=False):
    """
    Creates the OpenFOAM mesh from STL files.

    Args:
        case_dir:      Path to the OpenFOAM case directory
        stl_coldplate: Path to the full cold plate STL (Body_6)
        stl_fins:      Path to the fins STL (Body_3)
        pad:           Background mesh padding in meters
        nx, ny, nz:    Background mesh cell counts
        force:         If True, remesh even if mesh already exists
    """

    mesh_exists = os.path.exists(f"{case_dir}/constant/solid/polyMesh/boundary")
    if mesh_exists and not force:
        print("\nMesh already exists — skipping meshing. Pass force=True to remesh.")
        return

    print("\nRunning meshing pipeline...")
    os.makedirs(f"{case_dir}/constant/triSurface", exist_ok=True)
    os.makedirs(f"{case_dir}/constant/geometry", exist_ok=True)
    os.makedirs(f"{case_dir}/system", exist_ok=True)

    # ── Minimal controlDict required by all OpenFOAM utilities ────────────────
    write_file(f"{case_dir}/system/controlDict", """FoamFile
    {
        format      ascii;
        class       dictionary;
        object      controlDict;
    }
    application     foamMultiRun;
    startFrom       startTime;
    startTime       0;
    stopAt          endTime;
    endTime         1;
    deltaT          1;
    writeControl    timeStep;
    writeInterval   1;
    """)

    # ── Scale STL files to meters ──────────────────────────────────────────────
    for src, dst in [(stl_coldplate, "coldplate.stl"), (stl_fins, "fins.stl")]:
        mesh = pv.read(src)
        mesh.scale(0.001, inplace=True)
        mesh.save(f"{case_dir}/constant/triSurface/{dst}")
        mesh.save(f"{case_dir}/constant/geometry/{dst}")
        print(f"Saved {dst}, bounds: {[round(x,4) for x in mesh.bounds]}")

    # ── Compute bounds and locationInMesh ──────────────────────────────────────
    body6 = pv.read(f"{case_dir}/constant/triSurface/coldplate.stl")
    xmin, xmax = body6.bounds[0] - pad, body6.bounds[1] + pad
    ymin, ymax = body6.bounds[2] - pad, body6.bounds[3] + pad
    zmin, zmax = body6.bounds[4] - pad, body6.bounds[5] + pad
    cx, cy, cz = body6.center
    print(f"Bounds: X={xmin:.4f} to {xmax:.4f}, Y={ymin:.4f} to {ymax:.4f}, Z={zmin:.4f} to {zmax:.4f}")
    print(f"locationInMesh: ({cx:.6f} {cy:.6f} {cz:.6f})")

    # ── Write mesh dicts ───────────────────────────────────────────────────────
    write_file(f"{case_dir}/system/surfaceFeaturesDict", f"""FoamFile
{{
    format ascii; class dictionary; object surfaceFeaturesDict;
}}
surfaces ( "coldplate.stl" "fins.stl" );
includedAngle 150;
""")

    write_file(f"{case_dir}/system/blockMeshDict", f"""FoamFile
{{
    format ascii; class dictionary; object blockMeshDict;
}}
vertices
(
    ({xmin} {ymin} {zmin}) ({xmax} {ymin} {zmin})
    ({xmax} {ymax} {zmin}) ({xmin} {ymax} {zmin})
    ({xmin} {ymin} {zmax}) ({xmax} {ymin} {zmax})
    ({xmax} {ymax} {zmax}) ({xmin} {ymax} {zmax})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
defaultPatch {{ name default; type patch; }}
""")

    write_file(f"{case_dir}/system/snappyHexMeshDict", f"""FoamFile
{{
    format ascii; class dictionary; object snappyHexMeshDict;
}}
castellatedMesh true; snap true; addLayers false; mergeTolerance 1e-6;
geometry
{{
    coldplate {{ type triSurfaceMesh; file "coldplate.stl"; }}
    fins      {{ type triSurfaceMesh; file "fins.stl"; }}
}};
castellatedMeshControls
{{
    maxLocalCells 100000; maxGlobalCells 2000000;
    minRefinementCells 0; maxLoadUnbalance 0.10;
    nCellsBetweenLevels 2; resolveFeatureAngle 30;
    allowFreeStandingZoneFaces true;
    features
    (
        {{ file "coldplate.eMesh"; level 2; }}
        {{ file "fins.eMesh";      level 3; }}
    );
    refinementSurfaces
    {{
        coldplate {{ level (2 2); patchInfo {{ type wall; }} }}
        fins
        {{
            level (3 3); faceZone fins; cellZone solid;
            mode insidePoint;
            insidePoint ({cx:.6f} {cy/2:.6f} {cz:.6f});
        }}
    }}
    refinementRegions {{}}
    locationInMesh ({cx:.6f} {cy:.6f} {cz:.6f});
}}
snapControls {{ nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; }}
addLayersControls
{{
    relativeSizes true; layers {{}}
    expansionRatio 1.2; finalLayerThickness 0.3; minThickness 0.1;
}}
meshQualityControls
{{
    maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4;
    maxConcave 80; minFlatness 0.5; minVol 1e-13; minTetQuality -1e30;
    minArea -1; minTwist 0.02; minDeterminant 0.001; minFaceWeight 0.05;
    minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
}}
""")

    write_file(f"{case_dir}/system/solid/topoSetDict", """FoamFile
{
    format ascii; class dictionary; object topoSetDict;
}
actions
(
    { name inletFaces;  type faceSet; action new;    source patchToFace; patch default; }
    { name inletFaces;  type faceSet; action subset; source boxToFace;
      box (-0.046 -0.008 0.041) (0.046 0.036 0.083); }
    { name outletFaces; type faceSet; action new;    source patchToFace; patch default; }
    { name outletFaces; type faceSet; action subset; source boxToFace;
      box (-0.046 -0.008 -0.067) (0.046 0.036 -0.040); }
);
""")

    write_file(f"{case_dir}/system/solid/createPatchDict", """FoamFile
{
    format ascii; class dictionary; object createPatchDict;
}
pointSync false;
patches
(
    { name inlet;  patchInfo { type patch; } constructFrom set; set inletFaces;  }
    { name outlet; patchInfo { type patch; } constructFrom set; set outletFaces; }
);
""")

    # ── Run meshing commands ───────────────────────────────────────────────────
    run_cmd("source /opt/openfoam11/etc/bashrc && blockMesh", case_dir)
    run_cmd("source /opt/openfoam11/etc/bashrc && surfaceFeatures", case_dir)

    for emesh in glob.glob(f"{case_dir}/constant/geometry/*.eMesh"):
        shutil.copy(emesh, f"{case_dir}/constant/triSurface/")
        print(f"Copied {os.path.basename(emesh)} to triSurface")

    run_cmd("source /opt/openfoam11/etc/bashrc && snappyHexMesh -overwrite 2>&1 | tee snappy.log", case_dir)
    run_cmd("source /opt/openfoam11/etc/bashrc && splitMeshRegions -cellZones -overwrite", case_dir)
    run_cmd("source /opt/openfoam11/etc/bashrc && topoSet -region solid", case_dir)
    run_cmd("source /opt/openfoam11/etc/bashrc && createPatch -overwrite -region solid", case_dir)

    print("\nMeshing complete.")


def plot_convergence(case_dir,show_plots):
    """
    Parses solver.log and plots residuals for all solved fields.
    """
    log_path = f"{case_dir}/solver.log"
    if not os.path.exists(log_path):
        print("No solver.log found")
        return

    # ── Parse residuals from log ───────────────────────────────────────────────
    pattern = re.compile(
        r"smoothSolver|GAMG.*?Solving for (\w+).*?Initial residual = ([\d.e+-]+)",
    )
    # More robust pattern
    solve_pattern = re.compile(
        r"Solving for (\w+), Initial residual = ([\d.e+-]+), Final residual = ([\d.e+-]+)"
    )

    residuals = {}  # field -> list of initial residuals
    times     = {}  # field -> list of time steps

    current_time = None
    with open(log_path, 'r') as f:
        for line in f:
            # Detect time step
            time_match = re.match(r"^\s+Time = ([\d.]+)s?$", line)
            if time_match:
                current_time = float(time_match.group(1))
                continue

            # Detect residual lines
            solve_match = solve_pattern.search(line)
            if solve_match and current_time is not None:
                field   = solve_match.group(1)
                resid   = float(solve_match.group(2))
                if field not in residuals:
                    residuals[field] = []
                    times[field]     = []
                residuals[field].append(resid)
                times[field].append(current_time)

    if not residuals:
        print("No residuals found in log")
        return

    # ── Plot ──────────────────────────────────────────────────────────────────
    # Group fields for cleaner subplots
    fluid_fields = [f for f in residuals if f in ["Ux", "Uy", "Uz", "p_rgh", "h"]]
    solid_fields = [f for f in residuals if f in ["e"]]
    other_fields = [f for f in residuals if f not in fluid_fields + solid_fields]
    all_groups   = [("Fluid", fluid_fields), ("Solid", solid_fields)]
    if other_fields:
        all_groups.append(("Other", other_fields))

    fig, axes = plt.subplots(len(all_groups), 1, figsize=(12, 4 * len(all_groups)))
    if len(all_groups) == 1:
        axes = [axes]

    colors = plt.cm.tab10.colors

    for ax, (group_name, fields) in zip(axes, all_groups):
        for i, field in enumerate(fields):
            if field in residuals:
                ax.semilogy(
                    times[field], residuals[field],
                    label=field, color=colors[i % len(colors)], linewidth=1.5
                )
        ax.set_title(f"{group_name} Residuals")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Initial Residual")
        ax.legend(loc="upper right")
        ax.grid(True, which="both", alpha=0.3)
        ax.axhline(1e-4, color="gray", linestyle="--", alpha=0.5, label="1e-4 target")

    plt.suptitle("SIMPLE Convergence", fontsize=14)
    plt.tight_layout()

    plot_path = f"{case_dir}/convergence.png"
    plt.savefig(plot_path, dpi=150)
    print(f"Saved convergence plot: {plot_path}")
    if show_plots==True:
        plt.show()


# RUN SOLVER
def run_solver(case_dir, T_inlet, flow_rate, heat_flux, fin_volume, end_time=2500, write_interval=500,show_plots=True):
    """
    Writes all BC/property files and runs foamMultiRun.

    Args:
        case_dir:      Path to the OpenFOAM case directory
        T_inlet:       Inlet water temperature in Kelvin
        flow_rate:     Volumetric flow rate in m³/s
        heat_flux:     Total heat input in Watts
        fin_volume:    Fin volume in m³ (from checkMesh -region region1)
        end_time:      Number of SIMPLE iterations
        write_interval: How often to write results
    """

    pad        = 0.005
    A_end_face = (0.0401*2 + pad*2) * (0.0305 + 0.0023 + pad*2)
    U_eff      = flow_rate / A_end_face
    q_vol      = heat_flux / fin_volume
    A_base     = 0.080 * 0.138

    print(f"Inlet velocity: {U_eff:.6f} m/s")
    print(f"Volumetric heat source: {q_vol:.1f} W/m³")

    # controlDict 
    write_file(f"{case_dir}/system/controlDict", f"""FoamFile
{{
    format ascii; class dictionary; object controlDict;
}}
application foamMultiRun;
regionSolvers {{ solid fluid; region1 solid; }}
startFrom startTime; startTime 0;
stopAt endTime; endTime {end_time};
deltaT 1;
writeControl timeStep; writeInterval {write_interval};
purgeWrite 2; writeFormat binary; writePrecision 6;
writeCompression off; timeFormat general; timePrecision 6;
runTimeModifiable true;
functions
{{
    patchAverageInlet
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   {write_interval};
        operation       areaAverage;
        regionType      patch;
        name            inlet;
        fields          (p_rgh p);
        region          solid;
    }}
    patchAverageOutlet
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   {write_interval};
        operation       areaAverage;
        regionType      patch;
        name            outlet;
        fields          (p_rgh p);
        region          solid;
    }}
}}
""")

    # fvSchemes / fvSolution fluid
    write_file(f"{case_dir}/system/solid/fvSchemes", """FoamFile
{
    format ascii; class dictionary; object fvSchemes;
}
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      Gauss upwind;
    div(phi,h)      Gauss upwind;
    div(phi,K)      Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes    { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes       { default limited 0.5; }
""")

    write_file(f"{case_dir}/system/solid/fvSolution", """FoamFile
{
    format ascii; class dictionary; object fvSolution;
}
solvers
{
    p_rgh      { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    p_rghFinal { $p_rgh; relTol 0; }
    rho        { solver diagonal; }
    rhoFinal   { $rho; }
    "(U|h|k|epsilon)"      { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }
    "(U|h|k|epsilon)Final" { $U; relTol 0; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent yes;
    pRefCell 0; pRefValue 0;
}
relaxationFactors
{
    equations { U 0.3; h 0.5; }
    fields    { p_rgh 0.3; }
}
""")

    # fvSchemes / fvSolution solid 
    write_file(f"{case_dir}/system/region1/fvSchemes", """FoamFile
{
    format ascii; class dictionary; object fvSchemes;
}
ddtSchemes       { default steadyState; }
gradSchemes      { default Gauss linear; }
divSchemes       { default none; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes    { default limited 0.5; }
""")

    write_file(f"{case_dir}/system/region1/fvSolution", """FoamFile
{
    format ascii; class dictionary; object fvSolution;
}
solvers
{
    e      { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }
    eFinal { $e; relTol 0; }
}
SIMPLE { nNonOrthogonalCorrectors 2; }
relaxationFactors { equations { e 0.7; } }
""")

    # Physical properties 
    write_file(f"{case_dir}/constant/solid/physicalProperties", """FoamFile
{
    format ascii; class dictionary; object physicalProperties;
}
thermoType
{
    type heRhoThermo; mixture pureMixture; transport const;
    thermo hConst; equationOfState rhoConst; specie specie;
    energy sensibleEnthalpy;
}
mixture
{
    specie          { molWeight 18.0; }
    equationOfState { rho 998; }
    thermodynamics  { Cp 4182; Hf 0; }
    transport       { mu 1.002e-3; Pr 7.01; }
}
""")

    write_file(f"{case_dir}/constant/solid/momentumTransport", """FoamFile
{
    format ascii; class dictionary; object momentumTransport;
}
simulationType laminar;
""")

    write_file(f"{case_dir}/constant/solid/g", """FoamFile
{
    format ascii; class uniformDimensionedVectorField; object g;
}
dimensions [0 1 -2 0 0 0 0];
value (0 -9.81 0);
""")

    write_file(f"{case_dir}/constant/solid/pRef", """FoamFile
{
    format ascii; class uniformDimensionedScalarField; object pRef;
}
dimensions [1 -1 -2 0 0 0 0];
value 1e5;
""")

    write_file(f"{case_dir}/constant/region1/physicalProperties", """FoamFile
{
    format ascii; class dictionary; object physicalProperties;
}
thermoType
{
    type heSolidThermo; mixture pureMixture; transport constIsoSolid;
    thermo eConst; equationOfState rhoConst; specie specie;
    energy sensibleInternalEnergy;
}
mixture
{
    specie          { molWeight 63.5; }
    equationOfState { rho 8960; }
    thermodynamics  { Cv 385; Hf 0; }
    transport       { kappa 400; }
}
""")

    write_file(f"{case_dir}/constant/region1/fvModels", f"""FoamFile
{{
    format ascii; class dictionary; object fvModels;
}}
heatSource {{ type heatSource; select all; q {q_vol:.1f}; }}
""")

    write_file(f"{case_dir}/constant/region1/g", """FoamFile
{
    format ascii; class uniformDimensionedVectorField; object g;
}
dimensions [0 1 -2 0 0 0 0];
value (0 -9.81 0);
""")

    write_file(f"{case_dir}/system/fvSchemes", """FoamFile
{
    format ascii; class dictionary; object fvSchemes;
}
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
""")

    write_file(f"{case_dir}/system/fvSolution", """FoamFile
{
    format ascii; class dictionary; object fvSolution;
}
solvers {}
SIMPLE { nNonOrthogonalCorrectors 2; }
PIMPLE { nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 2; }
""")

    # Boundary conditions
    write_file(f"{case_dir}/0/solid/T", f"""FoamFile
{{
    format ascii; class volScalarField; object T;
}}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {T_inlet};
boundaryField
{{
    default          {{ type zeroGradient; }}
    inlet            {{ type fixedValue; value uniform {T_inlet}; }}
    outlet           {{ type zeroGradient; }}
    coldplate        {{ type zeroGradient; }}
    solid_to_region1 {{ type coupledTemperature; value uniform {T_inlet}; }}
}}
""")

    write_file(f"{case_dir}/0/solid/U", f"""FoamFile
{{
    format ascii; class volVectorField; object U;
}}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{{
    default          {{ type noSlip; }}
    inlet            {{ type fixedValue; value uniform (0 0 {-U_eff:.6f}); }}
    outlet           {{ type zeroGradient; }}
    coldplate        {{ type noSlip; }}
    solid_to_region1 {{ type noSlip; }}
}}
""")

    write_file(f"{case_dir}/0/solid/p_rgh", f"""FoamFile
{{
    format ascii; class volScalarField; object p_rgh;
}}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    default          {{ type fixedFluxPressure; value uniform 0; }}
    inlet            {{ type fixedFluxPressure; value uniform 0; }}
    outlet           {{ type fixedValue;        value uniform 0; }}
    coldplate        {{ type fixedFluxPressure; value uniform 0; }}
    solid_to_region1 {{ type fixedFluxPressure; value uniform 0; }}
}}
""")

    write_file(f"{case_dir}/0/solid/p", f"""FoamFile
{{
    format ascii; class volScalarField; object p;
}}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{{
    default          {{ type calculated; value uniform 1e5; }}
    inlet            {{ type calculated; value uniform 1e5; }}
    outlet           {{ type calculated; value uniform 1e5; }}
    coldplate        {{ type calculated; value uniform 1e5; }}
    solid_to_region1 {{ type calculated; value uniform 1e5; }}
}}
""")

    write_file(f"{case_dir}/0/region1/T", f"""FoamFile
{{
    format ascii; class volScalarField; object T;
}}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {T_inlet};
boundaryField
{{
    coldplate        {{ type zeroGradient; }}
    region1_to_solid {{ type coupledTemperature; value uniform {T_inlet}; }}
}}
""")

    # Clean old time dirs and run
    for time_dir in glob.glob(f"{case_dir}/[0-9]*/"):
        if not time_dir.endswith("/0/"):
            shutil.rmtree(time_dir)
            print(f"Removed: {time_dir}")

    for artefact in glob.glob(f"{case_dir}/0/solid/cellToRegion"):
        os.remove(artefact)
        print(f"Removed artefact: {artefact}")
    for artefact in glob.glob(f"{case_dir}/0/region1/cellToRegion"):
        os.remove(artefact)
        print(f"Removed artefact: {artefact}")

    run_cmd("source /opt/openfoam11/etc/bashrc && foamMultiRun 2>&1 | tee solver.log", case_dir)
    print("\nSolver complete.")
    plot_convergence(case_dir,show_plots)


# FUNCTION 3: VIEW RESULTS
def view_results(case_dir, T_inlet, heat_flux, show_plots=True):
    """
    Extracts thermal resistance, pressure drop, and optionally plots temperature fields.

    Args:
        case_dir:    Path to the OpenFOAM case directory
        T_inlet:     Inlet temperature in Kelvin
        heat_flux:   Total heat input in Watts
        show_plots:  If True, opens interactive pyvista windows

    Returns:
        dict with T_max, R_th, delta_p, and latest time step
    """

    # Find latest time step
    time_dirs = sorted([
        d for d in os.listdir(case_dir)
        if os.path.isdir(f"{case_dir}/{d}") and d.replace('.','').isdigit()
    ], key=float)
    latest = time_dirs[-1]
    print(f"Reading results from time: {latest}")

    # Convert to VTK
    for region in ["solid", "region1"]:
        subprocess.run(
            f"source /opt/openfoam11/etc/bashrc && foamToVTK -region {region} -time {latest}",
            shell=True, cwd=case_dir, env=OPENFOAM_ENV, executable="/bin/bash"
        )

    # Thermal resistance from solid region 
    T_max = None
    R_th  = None
    solid_vtk = glob.glob(f"{case_dir}/VTK/region1/coldplate_case_{latest}.vtk")
    if solid_vtk:
        solid = pv.read(solid_vtk[0])
        if "T" in solid.array_names:
            T_max = solid["T"].max()
            R_th  = (T_max - T_inlet) / heat_flux
            print(f"\n Thermal Results ")
            print(f"T_max  = {T_max:.2f} K  ({T_max-273.15:.2f}°C)")
            print(f"T_inlet= {T_inlet:.2f} K  ({T_inlet-273.15:.2f}°C)")
            print(f"R_th   = {R_th:.6f} K/W")

    # Pressure drop from inlet and outlet patches 
    delta_p = None
    inlet_vtk  = glob.glob(f"{case_dir}/VTK/solid/inlet/inlet_{latest}.vtk")
    outlet_vtk = glob.glob(f"{case_dir}/VTK/solid/outlet/outlet_{latest}.vtk")

    if inlet_vtk and outlet_vtk:
        inlet_mesh  = pv.read(inlet_vtk[0])
        outlet_mesh = pv.read(outlet_vtk[0])

        # Use p_rgh if available, fall back to p
        p_field = "p_rgh" if "p_rgh" in inlet_mesh.array_names else "p"

        if p_field in inlet_mesh.array_names and p_field in outlet_mesh.array_names:
            p_inlet  = inlet_mesh[p_field].mean()
            p_outlet = outlet_mesh[p_field].mean()
            delta_p  = p_inlet - p_outlet
            print(f"\n── Pressure Results ────────────────────")
            print(f"P_inlet  = {p_inlet:.2f} Pa")
            print(f"P_outlet = {p_outlet:.2f} Pa")
            print(f"ΔP       = {delta_p:.2f} Pa")
        else:
            print(f"Pressure field not found. Available: {inlet_mesh.array_names}")
    else:
        print("Inlet or outlet VTK patches not found")

    # ── Plots ─────────────────────────────────────────────────────────────────
    if show_plots:
        # Cold plate surface temperature
        coldplate_vtk = glob.glob(f"{case_dir}/VTK/region1/coldplate/coldplate_{latest}.vtk")
        if coldplate_vtk:
            cp = pv.read(coldplate_vtk[0])
            if "T" in cp.array_names:
                # Project onto full body6 STL
                body6_path = glob.glob(f"{case_dir}/constant/triSurface/coldplate.stl")
                if body6_path:
                    body6 = pv.read(body6_path[0])
                    body6_with_T = body6.sample(cp)
                    pl = pv.Plotter()
                    pl.add_mesh(body6_with_T, scalars="T", cmap="hot",
                                clim=[T_inlet, T_max] if T_max else None)
                    pl.add_scalar_bar("Temperature (K)")
                    pl.title = "Cold Plate Surface Temperature"
                    pl.show()

        # Fluid-solid interface temperature
        interface_vtk = glob.glob(f"{case_dir}/VTK/solid/solid_to_region1/solid_to_region1_{latest}.vtk")
        if interface_vtk:
            interface = pv.read(interface_vtk[0])
            if "T" in interface.array_names:
                pl = pv.Plotter()
                pl.add_mesh(interface, scalars="T", cmap="coolwarm")
                pl.add_scalar_bar("Temperature (K)")
                pl.title = "Fluid-Solid Interface Temperature"
                pl.show()

        print("\nAvailable VTK files in solid:")
        for root, dirs, files in os.walk(f"{case_dir}/VTK/solid"):
            for f in files:
                print(f"  {os.path.join(root, f)}")
    return {
        "T_max"   : T_max,
        "R_th"    : R_th,
        "delta_p" : delta_p,
        "time"    : latest
    }