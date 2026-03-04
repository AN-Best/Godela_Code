import os
import pyvista as pv
import subprocess
import shutil
import glob
import math

case_dir = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_case"
stl_folder = "/home/aaronbest/Desktop/GodelaChallenge/coldplate_inst_stl/coldplate_inst"

# ── Physical parameters ────────────────────────────────────────────────────────
T_inlet   = 313.15
flow_rate = 1.65e-3/60
heat_flux = 1100
pad       = 0.005

A_base      = 0.080 * 0.138
q_flux      = heat_flux / A_base
A_end_face  = (0.0401*2 + pad*2) * (0.0305 + 0.0023 + pad*2)
U_effective = flow_rate / A_end_face
fin_volume  = 4.23564e-06  # m³ from checkMesh
q_volumetric = heat_flux / fin_volume

print(f"Effective inlet velocity: {U_effective:.6f} m/s")
print(f"Heat flux: {q_flux:.1f} W/m²")
print(f"Volumetric heat source: {q_volumetric:.1f} W/m³")

# ── Create directories ─────────────────────────────────────────────────────────
os.makedirs(f"{case_dir}/constant/triSurface", exist_ok=True)
os.makedirs(f"{case_dir}/constant/geometry", exist_ok=True)

# ── Scale STL files to meters ──────────────────────────────────────────────────
for src, dst in [("Body_6.stl", "coldplate.stl"), ("Body_3.stl", "fins.stl")]:
    mesh = pv.read(f"{stl_folder}/{src}")
    mesh.scale(0.001, inplace=True)
    mesh.save(f"{case_dir}/constant/triSurface/{dst}")
    mesh.save(f"{case_dir}/constant/geometry/{dst}")
    print(f"Saved {dst} in meters, bounds: {[round(x,4) for x in mesh.bounds]}")

# ── Get bounds and center ──────────────────────────────────────────────────────
body6 = pv.read(f"{case_dir}/constant/triSurface/coldplate.stl")

xmin, xmax = body6.bounds[0] - pad, body6.bounds[1] + pad
ymin, ymax = body6.bounds[2] - pad, body6.bounds[3] + pad
zmin, zmax = body6.bounds[4] - pad, body6.bounds[5] + pad

nx, ny, nz = 50, 50, 50
cx, cy, cz = body6.center
loc_x, loc_y, loc_z = cx, cy, cz

print(f"Bounds (m): X={xmin:.4f} to {xmax:.4f}, Y={ymin:.4f} to {ymax:.4f}, Z={zmin:.4f} to {zmax:.4f}")
print(f"locationInMesh: ({loc_x:.6f} {loc_y:.6f} {loc_z:.6f})")

# ── surfaceFeaturesDict ────────────────────────────────────────────────────────
surface_features = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      surfaceFeaturesDict;
}
surfaces
(
    "coldplate.stl"
    "fins.stl"
);
includedAngle   150;
"""

# ── blockMeshDict ──────────────────────────────────────────────────────────────
block_mesh = f"""FoamFile
{{
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
vertices
(
    ({xmin} {ymin} {zmin})
    ({xmax} {ymin} {zmin})
    ({xmax} {ymax} {zmin})
    ({xmin} {ymax} {zmin})
    ({xmin} {ymin} {zmax})
    ({xmax} {ymin} {zmax})
    ({xmax} {ymax} {zmax})
    ({xmin} {ymax} {zmax})
);
blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);
defaultPatch
{{
    name default;
    type patch;
}}
"""

# ── snappyHexMeshDict ──────────────────────────────────────────────────────────
snappy = f"""FoamFile
{{
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}
castellatedMesh true;
snap            true;
addLayers       false;
mergeTolerance  1e-6;
geometry
{{
    coldplate {{ type triSurfaceMesh; file "coldplate.stl"; }}
    fins      {{ type triSurfaceMesh; file "fins.stl"; }}
}};
castellatedMeshControls
{{
    maxLocalCells       100000;
    maxGlobalCells      2000000;
    minRefinementCells  0;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 2;
    resolveFeatureAngle 30;
    allowFreeStandingZoneFaces true;
    features
    (
        {{ file "coldplate.eMesh"; level 2; }}
        {{ file "fins.eMesh"; level 3; }}
    );
    refinementSurfaces
    {{
        coldplate {{ level (2 2); patchInfo {{ type wall; }} }}
        fins
        {{
            level (3 3);
            faceZone  fins;
            cellZone  solid;
            mode      insidePoint;
            insidePoint ({loc_x:.6f} {loc_y/2:.6f} {loc_z:.6f});
        }}
    }}
    refinementRegions {{}}
    locationInMesh ({loc_x:.6f} {loc_y:.6f} {loc_z:.6f});
}}
snapControls
{{
    nSmoothPatch 3;
    tolerance    2.0;
    nSolveIter   30;
    nRelaxIter   5;
}}
addLayersControls
{{
    relativeSizes       true;
    layers              {{}}
    expansionRatio      1.2;
    finalLayerThickness 0.3;
    minThickness        0.1;
}}
meshQualityControls
{{
    maxNonOrtho         65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave          80;
    minFlatness         0.5;
    minVol              1e-13;
    minTetQuality       -1e30;
    minArea             -1;
    minTwist            0.02;
    minDeterminant      0.001;
    minFaceWeight       0.05;
    minVolRatio         0.01;
    minTriangleTwist    -1;
    nSmoothScale        4;
    errorReduction      0.75;
}}
"""

# ── controlDict ────────────────────────────────────────────────────────────────
control_dict = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      controlDict;
}
application     foamMultiRun;
regionSolvers
{
    solid       fluid;
    region1     solid;
}
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         2500;
deltaT          1;
writeControl    timeStep;
writeInterval   500;
purgeWrite      2;
writeFormat     binary;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
"""

# ── fvSchemes fluid ────────────────────────────────────────────────────────────
fv_schemes_fluid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      fvSchemes;
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
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
"""

# ── fvSolution fluid ───────────────────────────────────────────────────────────
fv_solution_fluid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    p_rgh
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.01;
    }
    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }
    rho
    {
        solver          diagonal;
    }
    rhoFinal
    {
        $rho;
    }
    "(U|h|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }
    "(U|h|k|epsilon)Final"
    {
        $U;
        relTol          0;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent               yes;
    pRefCell                 0;
    pRefValue                0;
}
relaxationFactors
{
    equations
    {
        U     0.3;
        h     0.5;
    }
    fields
    {
        p_rgh 0.3;
    }
}
"""

# ── fvSchemes solid ────────────────────────────────────────────────────────────
fv_schemes_solid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
"""

# ── fvSolution solid ───────────────────────────────────────────────────────────
fv_solution_solid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    e
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
    }
    eFinal
    {
        $e;
        relTol          0;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
}
relaxationFactors
{
    equations
    {
        e   0.7;
    }
}
"""

# ── physicalProperties fluid (water) ──────────────────────────────────────────
physical_properties_fluid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      physicalProperties;
}
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie          { molWeight 18.0; }
    equationOfState { rho 998; }
    thermodynamics  { Cp 4182; Hf 0; }
    transport       { mu 1.002e-3; Pr 7.01; }
}
"""

# ── physicalProperties solid (copper) ─────────────────────────────────────────
physical_properties_solid = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      physicalProperties;
}
thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIsoSolid;
    thermo          eConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}
mixture
{
    specie          { molWeight 63.5; }
    equationOfState { rho 8960; }
    thermodynamics  { Cv 385; Hf 0; }
    transport       { kappa 400; }
}
"""

# ── fvModels solid — volumetric heat source ────────────────────────────────────
fv_models_solid = f"""FoamFile
{{
    format      ascii;
    class       dictionary;
    object      fvModels;
}}
heatSource
{{
    type        heatSource;
    select      all;
    q           {q_volumetric:.1f};
}}
"""

# ── momentumTransport fluid ────────────────────────────────────────────────────
momentum_transport = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      momentumTransport;
}
simulationType laminar;
"""

# ── g ──────────────────────────────────────────────────────────────────────────
g_file = """FoamFile
{
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}
dimensions      [0 1 -2 0 0 0 0];
value           (0 -9.81 0);
"""

# ── pRef ──────────────────────────────────────────────────────────────────────
p_ref = """FoamFile
{
    format      ascii;
    class       uniformDimensionedScalarField;
    object      pRef;
}
dimensions      [1 -1 -2 0 0 0 0];
value           1e5;
"""

# ── BC: T fluid ────────────────────────────────────────────────────────────────
T_fluid = f"""FoamFile
{{
    format      ascii;
    class       volScalarField;
    object      T;
}}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {T_inlet};
boundaryField
{{
    default
    {{
        type            zeroGradient;
    }}
    inlet
    {{
        type            fixedValue;
        value           uniform {T_inlet};
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    coldplate
    {{
        type            zeroGradient;
    }}
    solid_to_region1
    {{
        type            coupledTemperature;
        value           uniform {T_inlet};
    }}
}}
"""

# ── BC: U fluid ────────────────────────────────────────────────────────────────
U_fluid = f"""FoamFile
{{
    format      ascii;
    class       volVectorField;
    object      U;
}}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{{
    default
    {{
        type            noSlip;
    }}
    inlet
    {{
        type            fixedValue;
        value           uniform (0 0 {-U_effective:.6f});
    }}
    outlet
    {{
        type            zeroGradient;
    }}
    coldplate
    {{
        type            noSlip;
    }}
    solid_to_region1
    {{
        type            noSlip;
    }}
}}
"""

# ── BC: p_rgh fluid ────────────────────────────────────────────────────────────
p_rgh_fluid = f"""FoamFile
{{
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    default         {{ type fixedFluxPressure; value uniform 0; }}
    inlet           {{ type fixedFluxPressure; value uniform 0; }}
    outlet          {{ type fixedValue;        value uniform 0; }}
    coldplate       {{ type fixedFluxPressure; value uniform 0; }}
    solid_to_region1 {{ type fixedFluxPressure; value uniform 0; }}
}}
"""

# ── BC: p fluid ────────────────────────────────────────────────────────────────
p_fluid = f"""FoamFile
{{
    format      ascii;
    class       volScalarField;
    object      p;
}}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{{
    default         {{ type calculated; value uniform 1e5; }}
    inlet           {{ type calculated; value uniform 1e5; }}
    outlet          {{ type calculated; value uniform 1e5; }}
    coldplate       {{ type calculated; value uniform 1e5; }}
    solid_to_region1 {{ type calculated; value uniform 1e5; }}
}}
"""

# ── BC: T solid ────────────────────────────────────────────────────────────────
T_solid = f"""FoamFile
{{
    format      ascii;
    class       volScalarField;
    object      T;
}}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {T_inlet};
boundaryField
{{
    coldplate
    {{
        type            zeroGradient;
    }}
    region1_to_solid
    {{
        type            coupledTemperature;
        value           uniform {T_inlet};
    }}
}}
"""

# ── topoSetDict ────────────────────────────────────────────────────────────────
topo_set_dict = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}
actions
(
    {
        name    inletFaces;
        type    faceSet;
        action  new;
        source  patchToFace;
        patch   default;
    }
    {
        name    inletFaces;
        type    faceSet;
        action  subset;
        source  boxToFace;
        box     (-0.046 -0.008 0.041) (0.046 0.036 0.083);
    }
    {
        name    outletFaces;
        type    faceSet;
        action  new;
        source  patchToFace;
        patch   default;
    }
    {
        name    outletFaces;
        type    faceSet;
        action  subset;
        source  boxToFace;
        box     (-0.046 -0.008 -0.067) (0.046 0.036 -0.040);
    }
);
"""

# ── createPatchDict ────────────────────────────────────────────────────────────
create_patch_dict = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      createPatchDict;
}
pointSync false;
patches
(
    {
        name inlet;
        patchInfo { type patch; }
        constructFrom set;
        set inletFaces;
    }
    {
        name outlet;
        patchInfo { type patch; }
        constructFrom set;
        set outletFaces;
    }
);
"""

# ── Write all files ────────────────────────────────────────────────────────────
def write_file(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    print(f"Written: {path}")

write_file(f"{case_dir}/system/surfaceFeaturesDict",                  surface_features)
write_file(f"{case_dir}/system/blockMeshDict",                        block_mesh)
write_file(f"{case_dir}/system/snappyHexMeshDict",                    snappy)
write_file(f"{case_dir}/system/controlDict",                          control_dict)
write_file(f"{case_dir}/system/solid/fvSchemes",                      fv_schemes_fluid)
write_file(f"{case_dir}/system/solid/fvSolution",                     fv_solution_fluid)
write_file(f"{case_dir}/system/solid/topoSetDict",                    topo_set_dict)
write_file(f"{case_dir}/system/solid/createPatchDict",                create_patch_dict)
write_file(f"{case_dir}/system/region1/fvSchemes",                    fv_schemes_solid)
write_file(f"{case_dir}/system/region1/fvSolution",                   fv_solution_solid)
write_file(f"{case_dir}/constant/solid/physicalProperties",           physical_properties_fluid)
write_file(f"{case_dir}/constant/solid/momentumTransport",            momentum_transport)
write_file(f"{case_dir}/constant/solid/g",                            g_file)
write_file(f"{case_dir}/constant/solid/pRef",                         p_ref)
write_file(f"{case_dir}/constant/region1/physicalProperties",         physical_properties_solid)
write_file(f"{case_dir}/constant/region1/fvModels",                   fv_models_solid)
write_file(f"{case_dir}/constant/region1/g",                          g_file)
write_file(f"{case_dir}/0/solid/T",                                   T_fluid)
write_file(f"{case_dir}/0/solid/U",                                   U_fluid)
write_file(f"{case_dir}/0/solid/p_rgh",                               p_rgh_fluid)
write_file(f"{case_dir}/0/solid/p",                                   p_fluid)
write_file(f"{case_dir}/0/region1/T",                                 T_solid)

# ── OpenFOAM runner ────────────────────────────────────────────────────────────
openfoam_env = os.environ.copy()
openfoam_env["PATH"] = "/opt/openfoam11/bin:" + openfoam_env["PATH"]
openfoam_env["WM_PROJECT_DIR"] = "/opt/openfoam11"

def run_cmd(cmd):
    print(f"\nRunning: {cmd}")
    result = subprocess.run(
        cmd, shell=True, cwd=case_dir,
        env=openfoam_env,
        executable="/bin/bash"
    )
    if result.returncode != 0:
        print(f"ERROR: {cmd} failed with return code {result.returncode}")
        exit(1)
    print(f"Done: {cmd}")

# ── Check if mesh already exists ──────────────────────────────────────────────
mesh_exists = os.path.exists(f"{case_dir}/constant/solid/polyMesh/boundary")

if mesh_exists:
    print("\nMesh already exists — skipping meshing steps")
else:
    print("\nNo mesh found — running full meshing pipeline")
    run_cmd("source /opt/openfoam11/etc/bashrc && blockMesh")
    run_cmd("source /opt/openfoam11/etc/bashrc && surfaceFeatures")

    for emesh in glob.glob(f"{case_dir}/constant/geometry/*.eMesh"):
        shutil.copy(emesh, f"{case_dir}/constant/triSurface/")
        print(f"Copied {os.path.basename(emesh)} to triSurface")

    run_cmd("source /opt/openfoam11/etc/bashrc && snappyHexMesh -overwrite 2>&1 | tee snappy.log")
    run_cmd("source /opt/openfoam11/etc/bashrc && splitMeshRegions -cellZones -overwrite")
    run_cmd("source /opt/openfoam11/etc/bashrc && topoSet -region solid")
    run_cmd("source /opt/openfoam11/etc/bashrc && createPatch -overwrite -region solid")

# ── Clean old time directories ─────────────────────────────────────────────────
for time_dir in glob.glob(f"{case_dir}/[0-9]*/"):
    if not time_dir.endswith("/0/"):
        shutil.rmtree(time_dir)
        print(f"Removed: {time_dir}")

run_cmd("source /opt/openfoam11/etc/bashrc && foamMultiRun 2>&1 | tee solver.log")