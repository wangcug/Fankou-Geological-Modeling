"""Evaluate Monte Carlo geological realizations against private CSAMT observations.

This implements the workflow stated in the accompanying manuscript: assign
representative lithological resistivities, apply anisotropic Gaussian smoothing,
interpolate measured apparent resistivity with IDW, and rank Pearson correlations.
The smoothing scales and IDW settings are user-configurable modeling assumptions:
the manuscript does not provide numerical Gaussian sigma/IDW hyperparameters.
This proxy comparison is NOT a full electromagnetic forward simulation.
"""

import argparse
import csv
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

from project_paths import required_file, output_path

# Table 1: use the mean of the representative/dominant lithology for each unit.
# Units with multiple lithologies require an interpretive choice; update the
# mapping if a different dominant lithology is defined by the private records.
UNIT_RESISTIVITY_OHM_M = {
    'C2+3ht': 22003.0,
    'D3tc': 400.0,
    'D3tb': 2538.0,
    'D3ta': 9789.0,
    'D2db': 5521.0,
    'D2da': 9789.0,
    'D1-2gt': 5521.0,
}


def cell_center_coordinates(extent, resolution):
    """Return voxel centers in GemPy/VTK flat order (X fastest)."""
    nx, ny, nz = resolution
    xs = extent[0] + (np.arange(nx) + 0.5) * (extent[1] - extent[0]) / nx
    ys = extent[2] + (np.arange(ny) + 0.5) * (extent[3] - extent[2]) / ny
    zs = extent[4] + (np.arange(nz) + 0.5) * (extent[5] - extent[4]) / nz
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing='ij')
    return np.column_stack([xx.ravel(order='F'), yy.ravel(order='F'), zz.ravel(order='F')])


def idw_interpolate(observation_xyz, observation_values, target_xyz, neighbors=12, power=2.0):
    """k-nearest-neighbor inverse-distance weighting, exact at sample locations."""
    tree = cKDTree(observation_xyz)
    out = np.empty(len(target_xyz), dtype=float)
    k = min(neighbors, len(observation_values))
    for lo in range(0, len(out), 25000):
        hi = min(lo + 25000, len(out))
        distances, indices = tree.query(target_xyz[lo:hi], k=k)
        distances = np.atleast_2d(distances) if k > 1 else distances[:, None]
        indices = np.atleast_2d(indices) if k > 1 else indices[:, None]
        exact = distances == 0
        weights = np.zeros_like(distances)
        nonexact_rows = ~exact.any(axis=1)
        weights[nonexact_rows] = 1 / np.maximum(distances[nonexact_rows], 1e-12) ** power
        weights[~nonexact_rows] = exact[~nonexact_rows].astype(float)
        out[lo:hi] = (weights * observation_values[indices]).sum(axis=1) / weights.sum(axis=1)
    return out


def lithology_to_resistivity(lith_flat, id_to_name):
    names = {int(k): v for k, v in id_to_name.items()}
    result = np.empty(len(lith_flat), dtype=float)
    for lith_id in np.unique(lith_flat):
        unit_name = names.get(int(lith_id))
        if unit_name not in UNIT_RESISTIVITY_OHM_M:
            raise ValueError(f'No lithological resistivity for voxel ID {lith_id} (unit={unit_name!r}). '
                             'Check the saved ID mapping and extend UNIT_RESISTIVITY_OHM_M if necessary.')
        result[lith_flat == lith_id] = UNIT_RESISTIVITY_OHM_M[unit_name]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--simulations', default=output_path('mc_simulation_results.pkl'))
    parser.add_argument('--observations', default=None, help='CSAMT Excel file with X,Y,Z,Rho columns')
    parser.add_argument('--sigma-horizontal-m', type=float, default=250.0,
                        help='Gaussian standard deviation in X and Y, metres (assumption; calibrate locally)')
    parser.add_argument('--sigma-vertical-m', type=float, default=25.0,
                        help='Gaussian standard deviation in Z, metres (assumption; calibrate locally)')
    parser.add_argument('--idw-neighbors', type=int, default=12)
    parser.add_argument('--idw-power', type=float, default=2.0)
    args = parser.parse_args()
    if args.sigma_horizontal_m <= 0 or args.sigma_vertical_m <= 0 or args.idw_neighbors < 1 or args.idw_power <= 0:
        parser.error('All smoothing scales, IDW neighbors and IDW power must be positive.')

    with open(args.simulations, 'rb') as stream:
        stored = pickle.load(stream)  # Load only trusted, locally generated files.
    params = stored['params']
    resolution = [int(v) for v in params['resolution']]
    if resolution != [42, 42, 42]:
        raise ValueError(f'Expected 42x42x42 simulation data, received {resolution}.')
    extent = np.array(params['extent'], dtype=float)
    id_to_name = stored['id_to_name']
    results = stored['simulation_results']
    if len(results) != 100:
        raise ValueError(f'Expected 100 realizations, got {len(results)}.')

    obs_path = args.observations or required_file('measured_csamt_resistivity.xlsx')
    observed = pd.read_excel(obs_path)
    observed.columns = observed.columns.astype(str).str.strip()
    missing = {'X', 'Y', 'Z', 'Rho'} - set(observed.columns)
    if missing:
        raise ValueError(f'Missing observation columns: {sorted(missing)}')
    observed = observed[['X', 'Y', 'Z', 'Rho']].apply(pd.to_numeric, errors='coerce').dropna()
    observed = observed[np.isfinite(observed.to_numpy()).all(axis=1) & (observed.Rho > 0)]
    if len(observed) < 3:
        raise ValueError('At least three valid CSAMT samples are required.')
    coords = cell_center_coordinates(extent, resolution)
    measured = idw_interpolate(observed[['X', 'Y', 'Z']].to_numpy(), observed.Rho.to_numpy(),
                               coords, args.idw_neighbors, args.idw_power)
    # Flattened GemPy voxel ordering: first X, then Y, then Z.
    dx = (extent[1] - extent[0]) / resolution[0]
    dy = (extent[3] - extent[2]) / resolution[1]
    dz = (extent[5] - extent[4]) / resolution[2]
    sigma_zyx = (args.sigma_vertical_m / dz, args.sigma_horizontal_m / dy,
                 args.sigma_horizontal_m / dx)
    ranking = []
    for row in results:
        lith = np.asarray(row['lith_block']).ravel()
        if lith.size != np.prod(resolution):
            raise ValueError('Realization grid does not match stored resolution.')
        rho = lithology_to_resistivity(lith, id_to_name).reshape(resolution[::-1])
        proxy = gaussian_filter(rho, sigma=sigma_zyx, mode='nearest').ravel()
        valid = np.isfinite(proxy) & np.isfinite(measured)
        if valid.sum() < 3 or np.std(proxy[valid]) == 0 or np.std(measured[valid]) == 0:
            corr = float('nan')
        else:
            corr = float(pearsonr(proxy[valid], measured[valid]).statistic)
        ranking.append({'realization': int(row['sim_num']) + 1,
                        'seed': int(row['seed']), 'pearson_r': corr})
        print(f"Realization {ranking[-1]['realization']:03d}: Pearson R = {corr:.4f}")
    ranking.sort(key=lambda r: np.nan_to_num(r['pearson_r'], nan=-np.inf), reverse=True)
    if not np.isfinite(ranking[0]['pearson_r']):
        raise RuntimeError('No finite model-to-measurement correlations were obtained.')
    dest = output_path('csamt_model_ranking.csv')
    with open(dest, 'w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=['realization', 'seed', 'pearson_r'])
        writer.writeheader()
        writer.writerows(ranking)
    print(f"Preferred realization: {ranking[0]['realization']} (Pearson R={ranking[0]['pearson_r']:.4f})")
    print(f'Results saved to: {dest}')


if __name__ == '__main__':
    main()
