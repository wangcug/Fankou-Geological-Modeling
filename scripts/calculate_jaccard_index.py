"""Compare 100 Monte Carlo voxel models against their 42^3 reference model.

The original study code uses matching-voxel fraction, labeled "Jaccard" in
its output. Conventional multi-class Jaccard (intersection over union) is
reported separately to prevent ambiguity; the original metric is retained.
"""
import argparse
import csv
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from project_paths import output_path


def compute_metrics(reference, candidate):
    """Return the paper-code agreement and macro classwise IoU."""
    ref, sim = np.asarray(reference).ravel(), np.asarray(candidate).ravel()
    if ref.shape != sim.shape or ref.size == 0:
        raise ValueError("The reference and realization must have equal nonzero size.")
    agreement = float(np.mean(ref == sim))
    ious = []
    for label in np.union1d(ref, sim):
        intersection = np.count_nonzero((ref == label) & (sim == label))
        union = np.count_nonzero((ref == label) | (sim == label))
        if union:
            ious.append(intersection / union)
    return agreement, float(np.mean(ious))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulations", default=output_path("mc_simulation_results.pkl"))
    args = parser.parse_args()
    with open(args.simulations, "rb") as stream:
        saved = pickle.load(stream)  # Only load your own trusted local results.
    reference = saved["reference_lith_block"]
    rows = []
    for result in saved["simulation_results"]:
        agreement, macro_iou = compute_metrics(reference, result["lith_block"])
        rows.append({"realization": int(result["sim_num"]) + 1,
                     "voxel_agreement_original_jaccard": agreement,
                     "macro_jaccard_iou": macro_iou})
    if not rows:
        raise ValueError("No Monte Carlo results were found.")
    destination = output_path("jaccard_scores.csv")
    with open(destination, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    values = np.array([r["voxel_agreement_original_jaccard"] for r in rows])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(values, bins=15, color="steelblue", edgecolor="white")
    ax.set_xlabel("Voxel agreement (original code's Jaccard metric)")
    ax.set_ylabel("Number of realizations")
    ax.set_title("Agreement of stochastic models with the reference model")
    ax.text(0.98, 0.98, f"N = {len(values)}\nMean = {values.mean():.4f}\n"
            f"Median = {np.median(values):.4f}\nRange = {values.min():.4f}–{values.max():.4f}",
            transform=ax.transAxes, va="top", ha="right")
    fig.tight_layout()
    fig.savefig(output_path("jaccard_distribution.png"), dpi=200)
    plt.close(fig)
    print(f"Wrote {destination}; N={len(rows)}, mean voxel agreement={values.mean():.4f}")


if __name__ == "__main__":
    main()
