# Software environment and dependencies

The table below reproduces the package versions provided by the authors. These versions are **not claimed to be a solved or independently tested environment**. A working Python environment with GemPy 2024.2.0 must be verified locally. Install using `pip install -r requirements.txt`, and if pip reports dependency conflicts, resolve them and document the working environment instead of silently changing package versions.

| Package | Version | Role in this study |
| --- | --- | --- |
| `numpy` | `1.24.3` | Arrays, voxel labels, numerical operations |
| `pandas` | `2.0.1` | Tables, geological coordinate observations and CSAMT input |
| `matplotlib` | `3.7.2` | 2D plots, histograms, cross-sections |
| `scipy` | `1.13.1` | KD-tree IDW, Gaussian filtering, Pearson correlation |
| `openpyxl` | `3.1.5` | Read private Excel workbooks through pandas |
| `geopandas` | `1.1.2` | Read GIS vectors and spatial coordinate layers |
| `fiona` | `1.10.1` | GIS file input/output backend |
| `shapely` | `2.1.1` | GIS geometry representation and iteration |
| `gempy` | `2024.2.0` | Implicit geological modeling and structural series |
| `gempy-viewer` | `2024.2.0` | GemPy interactive visualization support (not central to batch scripts) |
| `pyvista` | `0.45.0` | 3D voxel model rendering |
| `vtk` | `9.4.2` | Rendering engine used by PyVista |
| `jupyter` | `1.1.1` | Run the accompanying exploratory notebook |
| `nbformat` | `5.10.4` | Notebook file format and validation |
| `xarray` | `2024.7.0` | Scientific labeled arrays (environment dependency; not required by every script) |

## Suggested setup

```bash
python -m venv .venv
# Activate the environment (Windows: .venv\Scripts\activate; macOS/Linux: source .venv/bin/activate)
python -m pip install -r requirements.txt
python -m pip check
```

**Note:** This is a source release without the confidential source observations. Runtime validation of GemPy compatibility, geological inputs, and the final paper results is still required.
