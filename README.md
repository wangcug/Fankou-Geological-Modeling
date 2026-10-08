# Knowledge Graph-Guided 3D Geological Modeling Integrated with Geophysical Data

Research code for *Knowledge Graph-Guided 3D Geological Modeling Integrated with Geophysical Data: A Case Study of the Fankou Lead–Zinc Deposit, South China*.

## 1. System requirements and installation

- **Operating system:** Windows, Linux, or macOS with a working scientific Python environment. Three-dimensional visualization may require OpenGL/VTK graphics support; the rendering configuration is system-dependent.
- **Python:** Use a Python version compatible with the pinned scientific packages and GemPy 2024.2.0. The original Python interpreter version was not supplied, so no specific Python version is claimed to have been tested.
- **Dependencies:** All 15 package versions are listed in [`requirements.txt`](requirements.txt). Their purpose is explained in [`DEPENDENCIES.md`](DEPENDENCIES.md). These version pins are the author's recorded versions, not an independently verified clean installation.

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip check
```

## 2. Repository files and their functions

| Order | Code file | Purpose | Main input | Main output |
| --- | --- | --- | --- | --- |
| 1 | [`scripts/extract_geological_hard_data.py`](scripts/extract_geological_hard_data.py) | Extract coordinates from geological GIS features and assign DEM elevations for hard-data preparation | DEM; fault traces; optionally stratigraphic units and boundaries | Fault/boundary coordinate CSV files |
| 2 | [`scripts/build_3d_stratigraphic_fault_models.py`](scripts/build_3d_stratigraphic_fault_models.py) | Construct and visualize deterministic stratigraphic and fault models at 65 × 65 × 65 resolution | Reviewed `points.csv`, `orientations.csv` | 3D renderings and X–Z sections |
| 3 | [`3D_geological_model/3D_geological_model.ipynb`](3D_geological_model/3D_geological_model.ipynb) | Construct, inspect, and visualize the 3D geological model interactively | Reviewed geological modeling inputs | Interactive figures and notebook-generated results |
| 4 | [`scripts/monte_carlo_simulation.py`](scripts/monte_carlo_simulation.py) | Generate the 42 × 42 × 42 reference model and 100 constrained stochastic realizations | Reviewed `points.csv`, `orientations.csv` | `mc_simulation_results.pkl` |
| 5 | [`scripts/calculate_jaccard_index.py`](scripts/calculate_jaccard_index.py) | Compare stochastic voxel models with the reference model | `mc_simulation_results.pkl` | `jaccard_scores.csv` and histogram |
| 6 | [`scripts/calculate_pearson_r.py`](scripts/calculate_pearson_r.py) | Compare model-derived resistivity proxies with CSAMT observations and calculate Pearson R for model selection | Simulation results; measured CSAMT resistivity spreadsheet | `csamt_model_ranking.csv` |

[`scripts/project_paths.py`](scripts/project_paths.py) handles local input and output paths. It is a utility, not an additional research stage.

## 3. Workflow and execution order

The workflow follows the manuscript's Methods, Sections 3.2–3.6:

1. **Geological prior knowledge (Section 3.2).** Geological knowledge-graph relations define stratigraphic order and fault–stratum relationships. The knowledge-graph construction code is not included in this repository; the published model's units and faults are specified in the modeling scripts.
2. **Hard geological data extraction (Section 3.3).** Run `extract_geological_hard_data.py` to obtain mapped feature coordinates and DEM-derived elevations. Interpret cross-section observations separately and prepare the reviewed GemPy interface-point and orientation tables. Extracted coordinates alone are **not** equivalent to complete modeling inputs.
3. **3D geological modeling (Section 3.4).** Run `build_3d_stratigraphic_fault_models.py` for the deterministic stratigraphic/fault model. `3D_geological_model.ipynb` is a separate interactive route for inspecting the geological model; it is not required by the batch workflow.
4. **Monte Carlo uncertainty analysis (Section 3.5).** Run `monte_carlo_simulation.py` to generate a reference and 100 realizations, all at 42 × 42 × 42 resolution. Perturb stratigraphic interface points in Z (σ = 7 m) and fault points in X (σ = 15 m), retaining orientation constraints.
5. **Jaccard similarity calculation (Section 3.5).** Run `calculate_jaccard_index.py` on the Monte Carlo pickle. The original study code calls its voxel-label agreement measure “Jaccard”; the script also reports conventional classwise intersection-over-union separately.
6. **CSAMT-based evaluation (Section 3.6).** Run `calculate_pearson_r.py`, which assigns representative lithological resistivity, applies 3D Gaussian smoothing, interpolates measured CSAMT apparent resistivity by IDW, and ranks realizations using Pearson correlation R.

```text
Knowledge-graph topology + geological observations
                        |
     GIS + DEM --> extract_geological_hard_data.py
                        |
  Reviewed interface points + orientation constraints
                        |
            +-----------+----------------------+
            |                                  |
build_3d_stratigraphic_fault_models.py    monte_carlo_simulation.py
            |                                  |
     Deterministic 3D model                 100 realizations
                                               |
                                    +----------+----------+
                                    |                     |
                          calculate_jaccard_index.py  calculate_pearson_r.py
                                    |                     |
                            Model similarity        CSAMT evaluation / R
```

## 4. Required private input files

**No geological or geophysical datasets are included.** Supply your own files **outside the repository** and set `FANKOU_DATA_DIR` to that directory. By default, scripts look in `~/fankou_private_inputs` (your home directory); the default can be overridden.

| File name in your private input directory | Expected content |
| --- | --- |
| `dem.asc` | ESRI ASCII digital elevation model |
| `faults.shp` plus associated shapefile files | Mapped fault lines |
| `stratigraphic_units.shp` plus associated files (optional) | Geological unit polygons |
| `geological_boundaries.shp` plus associated files (optional) | Geological boundary lines |
| `points.csv` | Reviewed GemPy-compatible interface points, including spatial coordinates and geological surface identifiers |
| `orientations.csv` | Reviewed GemPy-compatible orientation observations |
| `measured_csamt_resistivity.xlsx` | Measured CSAMT apparent resistivity, with numeric `X`, `Y`, `Z`, and `Rho` columns (`Rho` in Ω·m) |

The geological input coordinates must share a consistent projected coordinate system in metres. The GIS extraction script does not generate the complete final GemPy interface and orientation datasets automatically.

Set the location of confidential files before running scripts:

```powershell
# Windows PowerShell (example path only)
$env:FANKOU_DATA_DIR = "D:\my_private_geological_inputs"
```

```bash
# Linux/macOS (example path only)
export FANKOU_DATA_DIR="$HOME/my_private_geological_inputs"
```

By default, generated results are written to a newly created `outputs/` directory. To change it, set `FANKOU_OUTPUT_DIR` to your preferred location. Generated results and confidential inputs are excluded from Git by `.gitignore`.

## 5. How to run

Run the following commands from the repository root, in order. Skip Step 1 when reviewed `points.csv` and `orientations.csv` are already available.

```bash
# Step 1 — Extract hard geological observations from private GIS and DEM inputs
python scripts/extract_geological_hard_data.py

# Step 2 — Build and visualize the deterministic 3D geological model (65 × 65 × 65)
python scripts/build_3d_stratigraphic_fault_models.py

# Step 3 — Generate 100 stochastic models and a reference model (42 × 42 × 42)
python scripts/monte_carlo_simulation.py

# Step 4 — Calculate model similarity / Jaccard-related statistics
python scripts/calculate_jaccard_index.py

# Step 5 — Calculate Pearson R against private CSAMT observations
python scripts/calculate_pearson_r.py
```

To open the 3D modeling notebook:

```bash
jupyter notebook 3D_geological_model/3D_geological_model.ipynb
```

The notebook is for interactive modeling and visualization and contains its own experimental calculations. The scripted sequence above is the recommended batch workflow.

## 6. Modeling settings and output interpretation

| Parameter | Value |
| --- | --- |
| Deterministic model resolution | 65 × 65 × 65 |
| Monte Carlo/reference resolution | 42 × 42 × 42 |
| Number of Monte Carlo realizations | 100 |
| Stratigraphic interface Z perturbation | Gaussian standard deviation of 7 m |
| Fault interface X perturbation | Gaussian standard deviation of 15 m |
| Study extent (X, Y, Z, m) | 462950–464050; 2777200–2779000; −800–100 |
| Modeled faults | F5, F6, F7 |

**Similarity metric:** The original similarity calculation counts identical lithology labels across corresponding voxels and divides by the number of voxels. This is voxel-label agreement, not the standard Jaccard/IoU formula. The script retains that value for continuity and separately reports conventional macro IoU.

**Pearson R:** The geophysical evaluation compares a smoothed lithology-based resistivity proxy against measured CSAMT apparent resistivity, not a full electromagnetic forward model. The numerical IDW and smoothing settings are configurable because their exact values are not fully specified in the manuscript.

**Execution note:** Syntax and selected helper-function tests can be performed without the restricted datasets, but a complete scientific rerun requires the correct confidential geological and CSAMT inputs and a compatible local GemPy installation.
