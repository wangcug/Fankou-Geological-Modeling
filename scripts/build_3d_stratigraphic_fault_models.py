from project_paths import OUTPUT_DIR as OUTPUT_DIR_PATH, required_file
"""
Fankou Pb-Zn deposit 3D geological modeling -  

"""
import os
import gc
import time
import numpy as np
import pyvista as pv
import vtk
import matplotlib
matplotlib.use('Agg')           
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm, to_hex
import matplotlib.font_manager as fm

pv.set_plot_theme('document')
import gempy as gp


np.random.seed(1515)

OUTPUT_DIR = str(OUTPUT_DIR_PATH)
os.makedirs(OUTPUT_DIR, exist_ok=True)
STRAT_PATH         = os.path.join(OUTPUT_DIR, "stratigraphic_model.png")
FAULT_PATH         = os.path.join(OUTPUT_DIR, "structural_model.png")
STRAT_SECTION_PATH = os.path.join(OUTPUT_DIR, "stratigraphy_xz_section.png")
FAULT_SECTION_PATH = os.path.join(OUTPUT_DIR, "faults_xz_section.png")

RESOLUTION  = [65, 65, 65]
CHUNK_SIZE  = 10000
WINDOW_SIZE = [1600, 1200]

ORIENTATIONS_PATH = required_file("orientations.csv")
POINTS_PATH = required_file("points.csv")

SECTION_Y = 2778100   



print("=" * 60)
t0 = time.time()

geo_model = gp.create_geomodel(
    project_name='Fankou',
    extent=[462950, 464050, 2777200, 2779000, -800, 100],
    resolution=RESOLUTION,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=ORIENTATIONS_PATH,
        path_to_surface_points=POINTS_PATH,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Fault_Series": ('F5', 'F6', 'F7'),
        "Strat_Series": ('C2+3ht', 'D3tc', 'D3tb', 'D3ta', 'D2db', 'D2da'),
    }
)
gp.set_is_fault(frame=geo_model.structural_frame, fault_groups=['Fault_Series'])

geo_model.interpolation_options.kernel_options.range  = 17
geo_model.interpolation_options.kernel_options.nugget = 0.015

opts = geo_model.interpolation_options.evaluation_options
for attr in ['evaluation_chunk_size', 'chunk_size', 'eval_chunk_size']:
    if hasattr(opts, attr):
        setattr(opts, attr, CHUNK_SIZE)
        break
if hasattr(opts, 'compute_scalar_gradient'):
    opts.compute_scalar_gradient = False

print(f"Computing model ({np.prod(RESOLUTION):,} voxels)...")
gp.compute_model(geo_model)
gc.collect()
print(f"Model computed in {time.time() - t0:.1f} s")


def build_imagedata(geo_model):
    e  = geo_model.grid.regular_grid.extent
    res = geo_model.grid.regular_grid.resolution
    nx, ny, nz = int(res[0]), int(res[1]), int(res[2])
    return pv.ImageData(
        dimensions=(nx + 1, ny + 1, nz + 1),
        spacing=((e[1]-e[0])/nx, (e[3]-e[2])/ny, (e[5]-e[4])/nz),
        origin=(e[0], e[2], e[4]),
    ), nx, ny, nz


extent = geo_model.grid.regular_grid.extent
xmin, xmax, ymin, ymax, zmin, zmax = extent
dx, dy, dz = xmax-xmin, ymax-ymin, zmax-zmin
cx, cy, cz  = (xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2
diag = max(dx, dy, dz)




CAMERA_POSITION_3D = [
    (cx - diag*0.9, cy - diag*1.2, cz + diag*0.6),  
    (cx, cy, cz),                                    
    (0, 0, 1),                                       
]


elements = list(geo_model.structural_frame.structural_elements)
FAULT_NAMES = {'F5', 'F6', 'F7'}

id2color = {}
id2name  = {}
for i, el in enumerate(elements):
    eid = i + 1
    id2color[eid] = el.color
    id2name[eid]  = el.name


lith       = np.asarray(geo_model.solutions.raw_arrays.lith_block).astype(int)
unique_ids = sorted(int(v) for v in np.unique(lith))


strat_ids = [u for u in unique_ids if id2name.get(u, '') not in FAULT_NAMES]
fault_ids = [u for u in unique_ids if id2name.get(u, '') in FAULT_NAMES]

n_strat  = len(strat_ids)
btr_cmap = plt.cm.RdBu_r                        
for i, sid in enumerate(strat_ids):
    frac = i / max(n_strat - 1, 1)              
    id2color[sid] = to_hex(btr_cmap(frac))
for fid in fault_ids:                            
    id2color[fid] = '#666666'

cmap_colors = [id2color.get(u, '#999999') for u in unique_ids]

annotations = {}
for u in unique_ids:
    name = id2name.get(u, f'id {u}')
    if name not in FAULT_NAMES:
        annotations[float(u)] = name

print("\n===== ID ->  /    =====")
for eid, name in id2name.items():
    print(f"  ID {eid}: {name}  {id2color[eid]}")

fb        = np.asarray(geo_model.solutions.raw_arrays.fault_block).astype(int)
unique_fb = sorted(int(v) for v in np.unique(fb))
BLOCK_COLORS = ['#1a3b6e', '#3a6ea5', '#6baed6', '#bdd7e7']

print("\n[1/4]  Stratigraphic model (3D) ...")
print(f"  lith_block unique IDs: {unique_ids}")

sargs_strat = dict(
    title='Stratigraphic units', title_font_size=16,
    position_x=0.02, position_y=0.05,
    width=0.10, height=0.88, vertical=True,
    n_labels=0, label_font_size=16, shadow=False, color='black',
)

mesh_lith, nx, ny, nz = build_imagedata(geo_model)
mesh_lith.cell_data['lith'] = lith

p1 = pv.Plotter(window_size=WINDOW_SIZE, title='Stratigraphic model')
p1.set_background('white')
p1.add_mesh(
    mesh_lith, scalars='lith',
    cmap=cmap_colors, n_colors=len(unique_ids),
    clim=[min(unique_ids)-0.5, max(unique_ids)+0.5],
    annotations=annotations, scalar_bar_args=sargs_strat,
    show_edges=False,
    ambient=0.30, diffuse=0.85, specular=0.15,
)
p1.camera_position = CAMERA_POSITION_3D

p1.show(screenshot=STRAT_PATH, auto_close=True)
print(f"Saved: {STRAT_PATH}")
del p1, mesh_lith
gc.collect()

print("\n[2/4]  Structural model (3D) ...")
print(f"  fault_block unique: {unique_fb}  -> {len(unique_fb)}  Fault block")

mesh_fault, _, _, _ = build_imagedata(geo_model)
mesh_fault.cell_data['fault_block'] = fb

p2 = pv.Plotter(window_size=WINDOW_SIZE, title='Structural model')
p2.set_background('white')
p2.add_mesh(
    mesh_fault, scalars='fault_block',
    cmap=BLOCK_COLORS[:len(unique_fb)], n_colors=len(unique_fb),
    clim=[min(unique_fb)-0.5, max(unique_fb)+0.5],
    show_edges=False, show_scalar_bar=False,
    ambient=0.30, diffuse=0.85, specular=0.15,
)
p2.camera_position = CAMERA_POSITION_3D

p2.show(screenshot=FAULT_PATH, auto_close=True)
print(f"Saved: {FAULT_PATH}")
del p2, mesh_fault
gc.collect()


y_frac = (SECTION_Y - ymin) / (ymax - ymin)
j_idx  = int(round(y_frac * ny))
j_idx  = max(0, min(j_idx, ny - 1))
y_actual = ymin + (j_idx + 0.5) * (ymax - ymin) / ny
print(f"XZ section: index {j_idx}/{ny}, actual Y={y_actual:.0f} m")

lith_3d = lith.reshape(nz, ny, nx)      # lith_3d[iz, iy, ix]
fb_3d   = fb.reshape(nz, ny, nx)

slice_lith = lith_3d[:, j_idx, :]
slice_fb   = fb_3d[:,   j_idx, :]

x_edges = np.linspace(xmin, xmax, nx + 1)
z_edges = np.linspace(zmin, zmax, nz + 1)

def save_section_mpl(data_2d, unique_vals, color_list, id_to_name,
                     exclude_names, save_path, title, section_y_val,
                     show_legend=True, fault_name_set=None):
    """
      matplotlib   2D XZ section 。
    data_2d  : shape (nz, nx)，  ID
    """
    fault_name_set = fault_name_set or set()

    mpl_cmap = ListedColormap(color_list)
    bounds   = [v - 0.5 for v in unique_vals] + [unique_vals[-1] + 0.5]
    norm     = BoundaryNorm(bounds, mpl_cmap.N)

    fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
    fig.patch.set_facecolor('white')

    pm = ax.pcolormesh(x_edges, z_edges, data_2d,
                       cmap=mpl_cmap, norm=norm, shading='flat')

    ax.set_xlabel('X coordinate (m)', fontsize=13)
    ax.set_ylabel('Z elevation (m)', fontsize=13)
    ax.set_title(f'{title}  (section Y ≈ {section_y_val:.0f} m)', fontsize=15)
    ax.tick_params(labelsize=11)
    ax.set_aspect('auto')

    
    if show_legend:
        patches = []
        for uid, color in zip(unique_vals, color_list):
            name = id_to_name.get(uid, f'id {uid}')
            if name in exclude_names:
                continue
            patches.append(mpatches.Patch(color=color, label=name))

        if patches:
            leg = ax.legend(
                handles=patches,
                loc='upper left',
                bbox_to_anchor=(-0.18, 1.0),   
                borderaxespad=0,
                fontsize=11,
                title='Stratigraphic units',
                title_fontsize=12,
                frameon=True,
                edgecolor='#888888',
            )
            fig.add_artist(leg)

    plt.tight_layout(rect=[0.18, 0, 1, 1])    
    plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✓  : {save_path}")


print(f"\n[3/4]  Stratigraphic model 2D XZ section (Y={SECTION_Y}) ...")

save_section_mpl(
    data_2d      = slice_lith,
    unique_vals  = unique_ids,
    color_list   = cmap_colors,
    id_to_name   = id2name,
    exclude_names= FAULT_NAMES,        
    save_path    = STRAT_SECTION_PATH,
    title        = 'Stratigraphic model XZ section',
    section_y_val= y_actual,
    show_legend  = True,
    fault_name_set=FAULT_NAMES,
)

print(f"\n[4/4]  Structural model 2D XZ section (Y={SECTION_Y}) ...")

fault_block_names = {u: f'Fault block {i+1}' for i, u in enumerate(unique_fb)}
fault_block_colors = BLOCK_COLORS[:len(unique_fb)]

fig, ax = plt.subplots(figsize=(14, 8), dpi=150)
fig.patch.set_facecolor('white')

mpl_cmap_fb = ListedColormap(fault_block_colors)
bounds_fb   = [v - 0.5 for v in unique_fb] + [unique_fb[-1] + 0.5]
norm_fb     = BoundaryNorm(bounds_fb, mpl_cmap_fb.N)

ax.pcolormesh(x_edges, z_edges, slice_fb,
              cmap=mpl_cmap_fb, norm=norm_fb, shading='flat')

ax.set_xlabel('X coordinate (m)', fontsize=13)
ax.set_ylabel('Z elevation (m)', fontsize=13)
ax.set_title(f'Structural model XZ section  (section Y = {y_actual:.0f} m)', fontsize=15)
ax.tick_params(labelsize=11)
ax.set_aspect('auto')

fb_patches = [mpatches.Patch(color=c, label=fault_block_names[u])
              for u, c in zip(unique_fb, fault_block_colors)]
ax.legend(
    handles=fb_patches,
    loc='upper left',
    bbox_to_anchor=(-0.18, 1.0),
    borderaxespad=0,
    fontsize=11,
    title='Fault block',
    title_fontsize=12,
    frameon=True,
    edgecolor='#888888',
)

plt.tight_layout(rect=[0.18, 0, 1, 1])
plt.savefig(FAULT_SECTION_PATH, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"  ✓  : {FAULT_SECTION_PATH}")

gc.collect()


total = time.time() - t0
print("\n" + "=" * 60)
print(f"✓  ！  {total:.1f}s ({total/60:.1f} min)")
print(f"\n   : {OUTPUT_DIR}")
print(f"  📊   (PyVista 3D,  → ,  ):")
print(f"     {STRAT_PATH}")
print(f"     {FAULT_PATH}")
print(f"  📈 XZ section (matplotlib 2D, Y ≈ {y_actual:.0f}):")
print(f"     {STRAT_SECTION_PATH}")
print(f"     {FAULT_SECTION_PATH}")
print("=" * 60)
print(f"\n💡  :   50   SECTION_Y  ")
print(f"     Y  : [{ymin}, {ymax}],   ≈ {(ymin+ymax)//2}")
print("=" * 60)
