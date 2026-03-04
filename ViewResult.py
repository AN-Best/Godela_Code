import pyvista as pv
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

case_dir  = "/home/aaronbest/Desktop/GodelaChallenge/Case2"
T_inlet   = 313.15
heat_flux = 1100
latest    = "2500"

# ── Load pedestal ──────────────────────────────────────────────────────────────
pedestal_file = glob.glob(
    f"{case_dir}/VTK/solid/solid_to_region1/solid_to_region1_{latest}.vtk")[0]
pedestal     = pv.read(pedestal_file)
pedestal_pts = pedestal.cell_data_to_point_data()
pts          = pedestal_pts.points
T_ped        = pedestal_pts["T"]
x_ped        = pts[:, 0]
z_ped        = pts[:, 2]
T_ped_min    = T_ped.min()
T_ped_max    = T_ped.max()
R_th         = (T_ped_max - T_inlet) / heat_flux

print(f"T_max = {T_ped_max:.2f} K  ({T_ped_max-273.15:.2f} °C)")
print(f"R_th  = {R_th:.6f} K/W")

# ── Load cold plate surface ────────────────────────────────────────────────────
coldplate_file = glob.glob(
    f"{case_dir}/VTK/solid/coldplate/coldplate_{latest}.vtk")[0]
coldplate_mesh = pv.read(coldplate_file)

# ── PyVista plot — cold plate 3D ───────────────────────────────────────────────
pl = pv.Plotter(window_size=(800, 600))
pl.add_mesh(coldplate_mesh, scalars="T", cmap="hot",
            clim=[T_inlet, coldplate_mesh["T"].max()],
            show_edges=False)
pl.add_scalar_bar("Temperature (K)")
pl.add_title(f"Cold Plate Surface  T_max={coldplate_mesh['T'].max()-273.15:.1f} °C")
pl.camera_position = "iso"
pl.show(screenshot="coldplate_surface.png")


fluid_file = glob.glob(f"{case_dir}/VTK/region1/Case2_{latest}.vtk")[0]
fluid_mesh = pv.read(fluid_file)

if "p" in fluid_mesh.array_names:
    bounds   = fluid_mesh.bounds
    z_min    = bounds[4]
    z_max    = bounds[5]
    z_range  = z_max - z_min

    # Slice at inlet (z_max) and outlet (z_min)
    inlet_slice  = fluid_mesh.slice(normal="z", origin=(0, 0, z_max - z_range*0.01))
    outlet_slice = fluid_mesh.slice(normal="z", origin=(0, 0, z_min + z_range*0.01))

    p_inlet  = inlet_slice["p"].mean()
    p_outlet = outlet_slice["p"].mean()
    delta_p  = p_inlet - p_outlet

    print(f"P_inlet  = {p_inlet:.4f} Pa")
    print(f"P_outlet = {p_outlet:.4f} Pa")
    print(f"ΔP       = {delta_p:.4f} Pa")
else:
    print(f"Available fields: {fluid_mesh.array_names}")

# ── Matplotlib — pedestal 2D heatmap ──────────────────────────────────────────
triang = mtri.Triangulation(x_ped, z_ped)

fig, ax = plt.subplots(figsize=(7, 6))
fig.patch.set_facecolor("white")
tpc = ax.tripcolor(triang, T_ped, cmap="hot",
                   vmin=T_ped_min, vmax=T_ped_max, shading="gouraud")
plt.colorbar(tpc, ax=ax, label="Temperature (K)")
ax.set_xlabel("X (m)")
ax.set_ylabel("Z (m)")
ax.set_title(f"Pedestal Heat Map\n"
             f"ΔT={T_ped_max-T_ped_min:.2f} K  |  R_th={R_th:.4f} K/W")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("pedestal_heatmap.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved pedestal_heatmap.png and coldplate_surface.png")