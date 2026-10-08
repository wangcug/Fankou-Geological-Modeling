"""Extract GIS-derived interface-point observations for 3D geological modeling.

This is a preprocessing stage only. It does not infer fault dip/orientation or
construct the final GemPy points.csv and orientations.csv automatically;
those require cross-section interpretation and geological review.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from project_paths import DATA_DIR, output_path


def read_ascii_dem(path):
    """Read an ESRI ASCII grid and return elevation, origin and pixel size."""
    with open(path, encoding="utf-8-sig") as stream:
        header = dict(stream.readline().strip().lower().split() for _ in range(6))
    ncols, nrows = int(header["ncols"]), int(header["nrows"])
    size = float(header["cellsize"])
    if "xllcorner" in header:
        x_left = float(header["xllcorner"])
    else:
        x_left = float(header["xllcenter"]) - size / 2
    if "yllcorner" in header:
        y_bottom = float(header["yllcorner"])
    else:
        y_bottom = float(header["yllcenter"]) - size / 2
    data = np.loadtxt(path, skiprows=6)
    if data.shape != (nrows, ncols):
        raise ValueError(f"DEM shape {data.shape} does not match {nrows}x{ncols} header")
    nodata = float(header.get("nodata_value", "-9999"))
    data[data == nodata] = np.nan
    return data, x_left, y_bottom + nrows * size, size


def dem_height(x, y, dem):
    """Nearest-cell elevation; unknown/out-of-domain elevations remain NaN."""
    elevation, x_left, y_top, size = dem
    col = int(np.floor((x - x_left) / size))
    row = int(np.floor((y_top - y) / size))
    if 0 <= row < elevation.shape[0] and 0 <= col < elevation.shape[1]:
        return float(elevation[row, col])
    return float("nan")


def geometry_lines(geometry):
    if geometry is None or geometry.is_empty:
        return
    if geometry.geom_type in ("LineString", "LinearRing"):
        yield geometry
    elif geometry.geom_type == "Polygon":
        yield geometry.exterior
        for ring in geometry.interiors:
            yield ring
    elif hasattr(geometry, "geoms"):
        for child in geometry.geoms:
            yield from geometry_lines(child)
    elif geometry.geom_type == "Point":
        yield geometry


def extract_points(path, feature_type, dem, target_crs=None):
    """Extract trace vertices and sample DEM elevations in the same projected CRS."""
    frame = gpd.read_file(path)
    if frame.crs is None:
        raise ValueError(f"Missing CRS: {path}. Specify a CRS in the source GIS files.")
    if target_crs and frame.crs != target_crs:
        frame = frame.to_crs(target_crs)
    rows = []
    names = [c for c in ("code", "CODE", "name", "NAME", "formation") if c in frame.columns]
    for index, feature in frame.iterrows():
        label = str(feature[names[0]]) if names and pd.notna(feature[names[0]]) else ""
        for part_index, line in enumerate(geometry_lines(feature.geometry)):
            for point_index, position in enumerate(line.coords):
                x, y = float(position[0]), float(position[1])
                rows.append({"feature_id": str(index), "part_index": part_index,
                             "point_index": point_index, "X": x, "Y": y,
                             "Z": dem_height(x, y, dem),
                             "feature_type": feature_type, "geological_unit": label})
    return pd.DataFrame(rows, columns=["feature_id", "part_index", "point_index",
                                      "X", "Y", "Z", "feature_type", "geological_unit"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dem", default=str(DATA_DIR / "dem.asc"))
    parser.add_argument("--faults", default=str(DATA_DIR / "faults.shp"))
    parser.add_argument("--units", default=str(DATA_DIR / "stratigraphic_units.shp"))
    parser.add_argument("--boundaries", default=str(DATA_DIR / "geological_boundaries.shp"))
    args = parser.parse_args()
    for filename in (args.dem, args.faults):
        if not Path(filename).is_file():
            parser.error(f"Missing required private input: {filename}")
    dem = read_ascii_dem(args.dem)
    layers = [(args.faults, "fault", "fault_coordinates.csv")]
    if Path(args.units).is_file():
        layers.append((args.units, "stratigraphic_unit", "stratigraphic_unit_coordinates.csv"))
    if Path(args.boundaries).is_file():
        layers.append((args.boundaries, "geological_boundary", "boundary_coordinates.csv"))
    results = []
    crs = None
    for filename, category, destination in layers:
        geo = gpd.read_file(filename)
        if geo.crs is None:
            raise ValueError(f"Missing CRS in {filename}")
        if not geo.crs.is_projected:
            raise ValueError(f"Projected CRS in metres required: {filename} ({geo.crs})")
        if crs is None:
            crs = geo.crs
        result = extract_points(filename, category, dem, crs)
        result.to_csv(output_path(destination), index=False)
        print(f"{category}: {len(result)} vertices -> {destination}")
        results.append(result)
    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(output_path("all_geological_coordinates.csv"), index=False)
    print(f"Saved {len(combined)} records. Verify raster CRS and geological interpretations before modeling.")


if __name__ == "__main__":
    main()
