"""Generate 100 topologically constrained stochastic geological models.

Inputs: private interface points and orientations. Outputs: a pickle containing
42x42x42 lithology labels, reference labels, model parameters and seeds.
Run calculate_jaccard_index.py and calculate_pearson_r.py afterward.
"""
from project_paths import OUTPUT_DIR as OUTPUT_DIR_PATH, required_file

import os
import gc
import time
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')           
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm, to_hex
import pandas as pd
import gempy as gp



MASTER_SEED = 1515


OUTPUT_DIR = str(OUTPUT_DIR_PATH)
ORIENTATIONS_PATH = required_file("orientations.csv")
POINTS_PATH = required_file("points.csv")


RESULTS_PICKLE = os.path.join(OUTPUT_DIR, "mc_simulation_results.pkl")


RESOLUTION = [42, 42, 42]
CHUNK_SIZE = 10000              


NUM_SIMS         = 100
STD_DEV_Z_STRATA = 7.0        
STD_DEV_X_FAULT  = 15.0        

SECTION_Y = 2778100            # section


FAULT_NAMES = {'F5', 'F6', 'F7'}   

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("\n" + "=" * 65)
print("Step 1: Compute the 42-cubed reference geological model...")
t_global = time.time()

np.random.seed(MASTER_SEED)    

geo_model = gp.create_geomodel(
    project_name='Fankou',
    extent=[462950, 464050, 2777200, 2779000, -800, 100],
    resolution=RESOLUTION,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=ORIENTATIONS_PATH,
        path_to_surface_points=POINTS_PATH,
    ),
)

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Fault_Series": ('F5', 'F6', 'F7'),
        "Strat_Series": ('C2+3ht', 'D3tc', 'D3tb', 'D3ta', 'D2db', 'D2da'),
    },
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

gp.compute_model(geo_model)
gc.collect()
print(f"Reference model ready in {time.time() - t_global:.1f} seconds")

elements = list(geo_model.structural_frame.structural_elements)
id2color, id2name = {}, {}
for i, el in enumerate(elements):
    eid = i + 1
    id2color[eid] = el.color
    id2name[eid]  = el.name


lith_ref   = np.asarray(geo_model.solutions.raw_arrays.lith_block).astype(int)
unique_ids = sorted(int(v) for v in np.unique(lith_ref))
strat_ids  = [u for u in unique_ids if id2name.get(u, '') not in FAULT_NAMES]
fault_ids  = [u for u in unique_ids if id2name.get(u, '') in FAULT_NAMES]

btr_cmap = plt.cm.RdBu_r
n_strat  = len(strat_ids)
for i, sid in enumerate(strat_ids):
    frac = i / max(n_strat - 1, 1)
    id2color[sid] = to_hex(btr_cmap(frac))
for fid in fault_ids:
    id2color[fid] = '#666666'

cmap_colors = [id2color.get(u, '#999999') for u in unique_ids]

print("\nStructural element IDs and colors:")
for eid, name in id2name.items():
    print(f"  ID {eid}: {name:12s}  {id2color[eid]}")

element_id_to_name = {el.id: el.name for el in geo_model.structural_frame.structural_elements}

original_sp  = geo_model.surface_points_copy.df.copy()
original_ori = geo_model.orientations_copy.df.copy()


def _add_element_col(df, id_col='id'):
    """  id  ，  surface / element  。"""
    df2 = df.copy()
    df2['surface'] = df2[id_col].map(element_id_to_name)
    df2['surface'] = df2['surface'].fillna(df2[id_col].astype(str))
    df2['element'] = df2['surface']
    return df2


sp_base  = _add_element_col(original_sp)
ori_base = _add_element_col(original_ori)

all_elements = list(dict.fromkeys(sp_base['element'].dropna().tolist() +
                                  ori_base['element'].dropna().tolist()))
FAULT_ELEMENTS = [e for e in all_elements if e in FAULT_NAMES]
STRAT_ELEMENTS = [e for e in all_elements
                  if e not in FAULT_ELEMENTS and e not in ('unknown', 'None')]

print(f"Fault elements: {FAULT_ELEMENTS}")
print(f"Stratigraphic elements: {STRAT_ELEMENTS}")


mapping_object = {
    "Fault_Series": ('F5', 'F6', 'F7'),
    "Strat_Series": ('C2+3ht', 'D3tc', 'D3tb', 'D3ta', 'D2db', 'D2da'),
}

extent = geo_model.grid.regular_grid.extent
x_min, x_max = float(extent[0]), float(extent[1])
z_min, z_max = float(extent[4]), float(extent[5])

print("\n" + "=" * 65)

if os.path.exists(RESULTS_PICKLE):
    
    print(f" ， : {RESULTS_PICKLE}")
    with open(RESULTS_PICKLE, 'rb') as f:
        saved_data = pickle.load(f)

    expected = {'num_sims': NUM_SIMS, 'master_seed': MASTER_SEED,
                'std_dev_z_strata': STD_DEV_Z_STRATA,
                'std_dev_x_fault': STD_DEV_X_FAULT,
                'resolution': RESOLUTION, 'extent': list(extent)}
    if saved_data.get('params') != expected or len(saved_data.get('simulation_results', [])) != NUM_SIMS:
        raise RuntimeError('Simulation cache parameters do not match this release. Delete mc_simulation_results.pkl from FANKOU_OUTPUT_DIR and rerun.')
    simulation_results = saved_data['simulation_results']
    if not np.array_equal(saved_data.get('reference_lith_block'), lith_ref.flatten()):
        raise RuntimeError('The reference voxel array differs from the cached run; delete the cache.')
    if saved_data.get('id_to_name') != id2name:
        raise RuntimeError('Cached lithology IDs do not match the current input model. Remove the cache and rerun.')
    print(f"Loaded {len(simulation_results)} cached realizations")
    print(f"Cached parameters: {saved_data.get('params', {})}")

else:
    
    print(f"Step 4     (n = {NUM_SIMS})")
    print(f"  σ_Z(strata) = {STD_DEV_Z_STRATA} m,   σ_X(fault) = {STD_DEV_X_FAULT} m")
    print(f"Grid resolution: {RESOLUTION}")
    t_mc = time.time()

    simulation_results = []
    for sim in range(NUM_SIMS):
        sim_seed = MASTER_SEED + sim * 17
        rng = np.random.default_rng(sim_seed)

        
        sp_pert  = sp_base.copy()
        ori_pert = ori_base.copy()

        
        mask_strat_sp  = sp_pert['element'].isin(STRAT_ELEMENTS)
        sp_pert.loc[mask_strat_sp, 'Z'] += rng.normal(0, STD_DEV_Z_STRATA, mask_strat_sp.sum())
        sp_pert.loc[mask_strat_sp, 'Z'] = sp_pert.loc[mask_strat_sp, 'Z'].clip(z_min, z_max)

        
        mask_fault_sp  = sp_pert['element'].isin(FAULT_ELEMENTS)
        sp_pert.loc[mask_fault_sp, 'X'] += rng.normal(0, STD_DEV_X_FAULT, mask_fault_sp.sum())
        sp_pert.loc[mask_fault_sp, 'X'] = sp_pert.loc[mask_fault_sp, 'X'].clip(x_min, x_max)

        tmp_sp_path  = os.path.join(OUTPUT_DIR, f'_tmp_sp_{sim}.csv')
        tmp_ori_path = os.path.join(OUTPUT_DIR, f'_tmp_ori_{sim}.csv')
        sp_pert.drop(columns=['element'], errors='ignore').to_csv(tmp_sp_path, index=False)
        ori_pert.drop(columns=['element'], errors='ignore').to_csv(tmp_ori_path, index=False)

        try:
            sim_model = gp.create_geomodel(
                project_name=f'Fankou_sim_{sim}',
                extent=[462950, 464050, 2777200, 2779000, -800, 100],
                resolution=RESOLUTION,
                importer_helper=gp.data.ImporterHelper(
                    path_to_orientations=tmp_ori_path,
                    path_to_surface_points=tmp_sp_path,
                ),
            )
            gp.map_stack_to_surfaces(gempy_model=sim_model, mapping_object=mapping_object)
            gp.set_is_fault(frame=sim_model.structural_frame, fault_groups=['Fault_Series'])

            sim_model.interpolation_options.kernel_options.range  = 17
            sim_model.interpolation_options.kernel_options.nugget = 0.015

            opts2 = sim_model.interpolation_options.evaluation_options
            for attr in ['evaluation_chunk_size', 'chunk_size', 'eval_chunk_size']:
                if hasattr(opts2, attr):
                    setattr(opts2, attr, CHUNK_SIZE)
                    break
            if hasattr(opts2, 'compute_scalar_gradient'):
                opts2.compute_scalar_gradient = False

            gp.compute_model(sim_model)

            lith_sim = np.asarray(sim_model.solutions.raw_arrays.lith_block).astype(int)
            simulation_results.append({
                'sim_num':    sim,
                'seed':       sim_seed,
                'lith_block': lith_sim.flatten(),
            })

            print(f"  Simulation {sim+1}/{NUM_SIMS} completed ({time.time() - t_mc:.0f}s)")

            del sim_model
            gc.collect()

        except Exception as e:
            print(f"  Simulation {sim+1}/{NUM_SIMS} failed: {e}")
        finally:
            for f in (tmp_sp_path, tmp_ori_path):
                if os.path.exists(f):
                    os.remove(f)

    print(f"Completed {len(simulation_results)}/{NUM_SIMS} realizations in {time.time() - t_mc:.1f}s")
    if len(simulation_results) != NUM_SIMS:
        raise RuntimeError("Not all Monte Carlo realizations succeeded; fix the input/model issue before evaluation.")

    saved_data = {
        'simulation_results': simulation_results,
        'id_to_name': id2name,
        'reference_lith_block': lith_ref.flatten(),
        'params': {
            'num_sims':         NUM_SIMS,
            'master_seed':      MASTER_SEED,
            'std_dev_z_strata': STD_DEV_Z_STRATA,
            'std_dev_x_fault':  STD_DEV_X_FAULT,
            'resolution':       RESOLUTION,
            'extent':           list(extent),
        },
    }
    with open(RESULTS_PICKLE, 'wb') as f:
        pickle.dump(saved_data, f, protocol=4)
    print(f'Saved simulation results: {RESULTS_PICKLE}')


print(f'Monte Carlo simulation completed: {len(simulation_results)} realizations.')
